#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

my $PI = 3.14159265358979323846;

# compress 5D hash code to 30-bit integer (6 bits per dimension)
sub compress_hash_key {
    my ($code) = @_;
    my $key_val = 0;
    for (my $i = 0; $i < 5; $i++) {
        my $val = $code->[$i];
        $val = 0 if $val < 0;
        $val = 63 if $val > 63;
        $key_val |= ($val << (6 * $i));
    }
    return $key_val;
}

sub clamp_int {
    my ($val, $lo, $hi) = @_;
    my $q = int($val >= 0.0 ? ($val + 0.5) : ($val - 0.5));
    $q = $lo if $q < $lo;
    $q = $hi if $q > $hi;
    return $q;
}

sub asin_clamp {
    my ($x) = @_;
    $x = 1.0 if $x > 1.0;
    $x = -1.0 if $x < -1.0;
    return atan2($x, sqrt(1.0 - $x * $x));
}

sub acos_clamp {
    my ($x) = @_;
    $x = 1.0 if $x > 1.0;
    $x = -1.0 if $x < -1.0;
    return atan2(sqrt(1.0 - $x * $x), $x);
}

sub main {
    my @dirs = File::Spec->splitdir(File::Spec->rel2abs($0));
    pop @dirs;
    pop @dirs;
    my $base_dir = File::Spec->catdir(@dirs);

    my $hyg_path = File::Spec->catfile($base_dir, "data", "hyg", "hyg_v44.csv");
    my $bin_db_path = File::Spec->catfile($base_dir, "data", "starfix_db.bin");

    print "Reading HYG Catalog from $hyg_path...\n";
    open(my $fh, '<', $hyg_path) or die "Cannot open catalog: $!";

    my $header = <$fh>;
    chomp $header;
    $header =~ s/"//g;
    my @cols = split(',', $header);
    my %col_map = map { $cols[$_] => $_ } 0 .. $#cols;

    my @stars;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line eq '';
        $line =~ s/"//g;
        my @row = split(',', $line);
        next if $row[$col_map{id}] eq '0';

        my $mag = $row[$col_map{mag}];
        next if $mag eq '' || $mag > 6.0;

        my $hip = $row[$col_map{hip}];
        $hip = -1 if !defined($hip) || $hip eq '';
        my $pmra = $row[$col_map{pmra}];
        $pmra = 0.0 if !defined($pmra) || $pmra eq '';
        my $pmdec = $row[$col_map{pmdec}];
        $pmdec = 0.0 if !defined($pmdec) || $pmdec eq '';

        push @stars, {
            id    => $row[$col_map{id}],
            hip   => int($hip),
            ra    => $row[$col_map{rarad}] + 0.0,
            dec   => $row[$col_map{decrad}] + 0.0,
            mag   => $mag + 0.0,
            pmra  => $pmra + 0.0,
            pmdec => $pmdec + 0.0,
        };
    }
    close($fh);

    print "Loaded " . scalar(@stars) . " stars with magnitude <= 6.0\n";
    @stars = sort { $a->{mag} <=> $b->{mag} } @stars;

    my @vectors;
    foreach my $star (@stars) {
        my ($ra, $dec) = ($star->{ra}, $star->{dec});
        push @vectors, [ cos($dec) * cos($ra), cos($dec) * sin($ra), sin($dec) ];
    }

    # merge only unresolved visual binaries (< 0.05 deg = 3 arcmin) in the star catalog
    my $cat_spacing_deg = 0.05;
    my $cos_cat_spacing = cos($cat_spacing_deg * $PI / 180.0);
    my @catalog_stars;
    my @catalog_vectors;
    for (my $i = 0; $i < scalar(@stars); $i++) {
        my $vec = $vectors[$i];
        my $blended = 0;
        foreach my $sel (@catalog_vectors) {
            my $dot = $vec->[0] * $sel->[0] + $vec->[1] * $sel->[1] + $vec->[2] * $sel->[2];
            if ($dot >= $cos_cat_spacing) {
                $blended = 1;
                last;
            }
        }
        if (!$blended) {
            push @catalog_stars, $stars[$i];
            push @catalog_vectors, $vec;
        }
    }

    my $num_stars = scalar(@catalog_stars);
    print "Kept $num_stars catalog stars after $cat_spacing_deg deg blend filter.\n";

    # select well-spaced seed subset (>= 1.4 deg) solely for quad hash indexing
    my $seed_spacing_deg = 1.4;
    my $cos_seed_spacing = cos($seed_spacing_deg * $PI / 180.0);
    my @seed_cat_idx;
    my @seed_vectors;
    for (my $i = 0; $i < $num_stars; $i++) {
        my $vec = $catalog_vectors[$i];
        my $too_close = 0;
        foreach my $sel (@seed_vectors) {
            my $dot = $vec->[0] * $sel->[0] + $vec->[1] * $sel->[1] + $vec->[2] * $sel->[2];
            if ($dot >= $cos_seed_spacing) {
                $too_close = 1;
                last;
            }
        }
        if (!$too_close) {
            push @seed_cat_idx, $i;
            push @seed_vectors, $vec;
        }
    }
    my $num_seeds = scalar(@seed_cat_idx);
    print "Selected $num_seeds quad seed stars (spacing >= $seed_spacing_deg deg).\n";

    my $max_fov_deg = 22.0;
    my $cos_max_fov = cos($max_fov_deg * $PI / 180.0);
    my $bin_factor = 50;

    print "Generating unique 4-star quads (max FOV = $max_fov_deg degrees)...\n";
    my %unique_quads;
    for (my $i = 0; $i < $num_seeds; $i++) {
        my $v_i = $seed_vectors[$i];
        my @neighbors;
        for (my $j = $i + 1; $j < $num_seeds; $j++) {
            my $v_j = $seed_vectors[$j];
            my $dot = $v_i->[0] * $v_j->[0] + $v_i->[1] * $v_j->[1] + $v_i->[2] * $v_j->[2];
            if ($dot >= $cos_max_fov) {
                push @neighbors, { cat_idx => $seed_cat_idx[$j], dist => acos_clamp($dot) };
            }
        }
        @neighbors = sort { $a->{dist} <=> $b->{dist} } @neighbors;
        splice(@neighbors, 10) if scalar(@neighbors) > 10;

        my $n = scalar(@neighbors);
        if ($n >= 3) {
            my $ci = $seed_cat_idx[$i];
            for (my $n1 = 0; $n1 < $n; $n1++) {
                for (my $n2 = $n1 + 1; $n2 < $n; $n2++) {
                    for (my $n3 = $n2 + 1; $n3 < $n; $n3++) {
                        my @q = sort { $a <=> $b } (
                            $ci,
                            $neighbors[$n1]->{cat_idx},
                            $neighbors[$n2]->{cat_idx},
                            $neighbors[$n3]->{cat_idx}
                        );
                        $unique_quads{ join(",", @q) } = \@q;
                    }
                }
            }
        }
    }

    my @quads = values %unique_quads;
    my $total_quads = scalar(@quads);
    print "Found $total_quads unique 4-star quads.\n";

    print "Hashing quads...\n";
    my @entry_records;
    foreach my $q (@quads) {
        my @v = map { $catalog_vectors[$_] } @$q;
        my @angles;
        for (my $a = 0; $a < 4; $a++) {
            for (my $b = $a + 1; $b < 4; $b++) {
                my $dx = $v[$a]->[0] - $v[$b]->[0];
                my $dy = $v[$a]->[1] - $v[$b]->[1];
                my $dz = $v[$a]->[2] - $v[$b]->[2];
                my $dist = sqrt($dx * $dx + $dy * $dy + $dz * $dz);
                push @angles, 2.0 * asin_clamp(0.5 * $dist);
            }
        }
        @angles = sort { $a <=> $b } @angles;
        my $largest = $angles[5];
        next if $largest < 1e-5;

        my @hash_code = map { int(($angles[$_] / $largest) * $bin_factor) } 0 .. 4;
        push @entry_records, {
            hash_key => compress_hash_key(\@hash_code),
            stars    => $q,
        };
    }

    print "Sorting entries...\n";
    @entry_records = sort {
        $a->{hash_key} <=> $b->{hash_key}
            || $a->{stars}[0] <=> $b->{stars}[0]
            || $a->{stars}[1] <=> $b->{stars}[1]
            || $a->{stars}[2] <=> $b->{stars}[2]
            || $a->{stars}[3] <=> $b->{stars}[3]
    } @entry_records;

    print "Writing binary database to $bin_db_path...\n";
    open(my $ofh, '>:raw', $bin_db_path) or die "Cannot open output file: $!";
    print $ofh "SFIX";
    print $ofh pack("V V V", $num_stars, scalar(@entry_records), $bin_factor);

    # 18-byte star record: int32 hip, uint32 ra_q, uint32 dec_q, int16 mag_mmag, int16 pm_ra_q, int16 pm_dec_q
    my $u32_max = 4294967295.0;
    foreach my $star (@catalog_stars) {
        my $ra_val = $star->{ra};
        $ra_val += 2.0 * $PI while $ra_val < 0.0;
        $ra_val -= 2.0 * $PI while $ra_val >= 2.0 * $PI;
        my $ra_q = int($ra_val * $u32_max / (2.0 * $PI) + 0.5);
        $ra_q = 4294967295 if $ra_q > 4294967295;

        my $dec_val = $star->{dec};
        $dec_val = -$PI / 2.0 if $dec_val < -$PI / 2.0;
        $dec_val = $PI / 2.0 if $dec_val > $PI / 2.0;
        my $dec_q = int(($dec_val + $PI / 2.0) * $u32_max / $PI + 0.5);
        $dec_q = 4294967295 if $dec_q > 4294967295;

        my $mag_q = clamp_int($star->{mag} * 1000.0, -32768, 32767);
        my $pm_ra_q = clamp_int($star->{pmra} * 5.0, -32768, 32767);
        my $pm_dec_q = clamp_int($star->{pmdec} * 5.0, -32768, 32767);

        print $ofh pack("l< V V s< s< s<", $star->{hip}, $ra_q, $dec_q, $mag_q, $pm_ra_q, $pm_dec_q);
    }

    foreach my $entry (@entry_records) {
        print $ofh pack("V v v v v", $entry->{hash_key}, @{ $entry->{stars} });
    }

    close($ofh);
    printf("Binary database compiled successfully! Total size: %.1f KB\n", (-s $bin_db_path) / 1024.0);
}

main();
