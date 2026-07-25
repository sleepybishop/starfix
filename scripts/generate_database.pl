#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;

# compress 5D hash code to 30-bit integer (6 bits per dimension)
sub compress_hash_key {
    my $code = shift;
    my $key_val = 0;
    for (my $i = 0; $i < 5; $i++) {
        my $val = $code->[$i];
        $val = 0 if $val < 0;
        $val = 63 if $val > 63;
        $key_val += $val * (64 ** $i);
    }
    return $key_val;
}

sub main {
    # find project base dir
    my $base_dir = File::Spec->rel2abs(File::Spec->path());
    my @dirs = File::Spec->splitdir(File::Spec->rel2abs($0));
    pop @dirs; # pop filename
    pop @dirs; # pop scripts dir
    $base_dir = File::Spec->catdir(@dirs);

    my $hyg_path = File::Spec->catfile($base_dir, "data", "hyg", "hyg_v44.csv");
    my $bin_db_path = File::Spec->catfile($base_dir, "data", "starfix_db.bin");

    print "Reading HYG Catalog from $hyg_path...\n";
    open(my $fh, '<', $hyg_path) or die "Cannot open catalog: $!";
    
    my $header = <$fh>;
    chomp $header;
    $header =~ s/"//g;
    my @cols = split(',', $header);
    my %col_map = map { $cols[$_] => $_ } 0..$#cols;
    
    my @stars;
    while (my $line = <$fh>) {
        chomp $line;
        next if $line eq '';
        $line =~ s/"//g;
        my @row = split(',', $line);
        next if $row[$col_map{id}] eq '0'; # skip sun
        
        my $mag = $row[$col_map{mag}];
        next if $mag eq '' || $mag > 6.0;
        
        my $hip = $row[$col_map{hip}];
        $hip = -1 if $hip eq '';
        
        push @stars, {
            id => $row[$col_map{id}],
            hip => $hip,
            proper => $row[$col_map{proper}] || "",
            ra => $row[$col_map{rarad}],
            dec => $row[$col_map{decrad}],
            mag => $mag
        };
    }
    close($fh);
    
    print "Loaded " . scalar(@stars) . " stars with magnitude <= 6.0\n";
    
    # sort stars by magnitude (brightest first)
    @stars = sort { $a->{mag} <=> $b->{mag} } @stars;
    
    # compute 3D unit vectors
    my @vectors;
    foreach my $star (@stars) {
        my $ra = $star->{ra};
        my $dec = $star->{dec};
        push @vectors, [
            cos($dec) * cos($ra),
            cos($dec) * sin($ra),
            sin($dec)
        ];
    }
    
    # spacing filter: keep brightest, space out by >= 1.5 degrees
    my $min_spacing_deg = 1.5;
    my $min_spacing_rad = $min_spacing_deg * 3.1415926535897932 / 180.0;
    my $cos_spacing = cos($min_spacing_rad);
    
    print "Filtering stars with a minimum spacing of $min_spacing_deg degrees...\n";
    my @filtered_stars;
    my @star_vectors;
    
    for (my $i = 0; $i < scalar(@stars); $i++) {
        my $vec = $vectors[$i];
        my $too_close = 0;
        foreach my $sel_vec (@star_vectors) {
            my $dot = $vec->[0]*$sel_vec->[0] + $vec->[1]*$sel_vec->[1] + $vec->[2]*$sel_vec->[2];
            if ($dot >= $cos_spacing) {
                $too_close = 1;
                last;
            }
        }
        if (!$too_close) {
            push @filtered_stars, $stars[$i];
            push @star_vectors, $vec;
        }
    }
    
    my $num_stars = scalar(@filtered_stars);
    print "Kept $num_stars stars after spacing filter.\n";
    
    # generate star patterns (max FOV = 15.0 deg)
    my $max_fov_deg = 15.0;
    my $max_fov_rad = $max_fov_deg * 3.1415926535897932 / 180.0;
    my $cos_max_fov = cos($max_fov_rad);
    my $bin_factor = 50;
    
    print "Generating unique 4-star quads (max FOV = $max_fov_deg degrees)...\n";
    my %unique_quads;
    
    for (my $i = 0; $i < $num_stars; $i++) {
        my $v_i = $star_vectors[$i];
        
        # find neighbors
        my @neighbors;
        for (my $j = $i + 1; $j < $num_stars; $j++) {
            my $v_j = $star_vectors[$j];
            my $dot = $v_i->[0]*$v_j->[0] + $v_i->[1]*$v_j->[1] + $v_i->[2]*$v_j->[2];
            if ($dot >= $cos_max_fov) {
                my $dist = acos_clamp($dot);
                push @neighbors, { idx => $j, dist => $dist };
            }
        }
        
        # sort neighbors by distance, keep 10 closest to prevent hash table bloat
        @neighbors = sort { $a->{dist} <=> $b->{dist} } @neighbors;
        if (scalar(@neighbors) > 10) {
            splice(@neighbors, 10);
        }
        
        # generate combinations
        if (scalar(@neighbors) >= 3) {
            my $n = scalar(@neighbors);
            for (my $n1 = 0; $n1 < $n; $n1++) {
                for (my $n2 = $n1 + 1; $n2 < $n; $n2++) {
                    for (my $n3 = $n2 + 1; $n3 < $n; $n3++) {
                        my @q = sort { $a <=> $b } ($i, $neighbors[$n1]->{idx}, $neighbors[$n2]->{idx}, $neighbors[$n3]->{idx});
                        my $key = join(",", @q);
                        $unique_quads{$key} = \@q;
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
    my $quad_idx = 0;
    
    foreach my $q (@quads) {
        if ($quad_idx % 20000 == 0 && $quad_idx > 0) {
            print "  Processed $quad_idx/$total_quads quads...\n";
        }
        $quad_idx++;
        
        # extract unit vectors
        my @v = (
            $star_vectors[$q->[0]],
            $star_vectors[$q->[1]],
            $star_vectors[$q->[2]],
            $star_vectors[$q->[3]]
        );
        
        # 6 pairwise angles
        my @angles;
        for (my $a = 0; $a < 4; $a++) {
            for (my $b = $a + 1; $b < 4; $b++) {
                my $dx = $v[$a]->[0] - $v[$b]->[0];
                my $dy = $v[$a]->[1] - $v[$b]->[1];
                my $dz = $v[$a]->[2] - $v[$b]->[2];
                my $dist = sqrt($dx*$dx + $dy*$dy + $dz*$dz);
                my $val = 0.5 * $dist;
                $val = 1.0 if $val > 1.0;
                my $angle = 2.0 * asin($val);
                push @angles, $angle;
            }
        }
        
        @angles = sort { $a <=> $b } @angles;
        my $largest_angle = $angles[5];
        next if $largest_angle < 1e-5;
        
        # compute ratios and bin factor
        my @hash_code;
        for (my $k = 0; $k < 5; $k++) {
            my $ratio = $angles[$k] / $largest_angle;
            push @hash_code, int($ratio * $bin_factor);
        }
        
        my $hash_key = compress_hash_key(\@hash_code);
        push @entry_records, {
            hash_key => $hash_key,
            stars => $q
        };
    }
    
    # Sort entries by hash_key ascending (critical for binary search!)
    print "Sorting entries...\n";
    @entry_records = sort { $a->{hash_key} <=> $b->{hash_key} } @entry_records;
    
    # Write binary database
    print "Writing binary database to $bin_db_path...\n";
    open(my $ofh, '>:raw', $bin_db_path) or die "Cannot open output file: $!";
    
    # Header: Magic, num_stars, num_entries, bin_factor
    print $ofh "SFIX";
    print $ofh pack("V V V", $num_stars, scalar(@entry_records), $bin_factor);
    
    # write star catalog records (signed 32-bit hip, quantized 16-bit ra, dec, quantized 8-bit mag)
    foreach my $star (@filtered_stars) {
        my $ra_val = $star->{ra};
        while ($ra_val < 0.0) { $ra_val += 2.0 * 3.1415926535897932; }
        while ($ra_val >= 2.0 * 3.1415926535897932) { $ra_val -= 2.0 * 3.1415926535897932; }
        my $ra_q = int($ra_val * 65535.0 / (2.0 * 3.1415926535897932) + 0.5);

        my $dec_val = $star->{dec};
        if ($dec_val < -3.1415926535897932 / 2.0) { $dec_val = -3.1415926535897932 / 2.0; }
        if ($dec_val > 3.1415926535897932 / 2.0) { $dec_val = 3.1415926535897932 / 2.0; }
        my $dec_q = int(($dec_val + 3.1415926535897932 / 2.0) * 65535.0 / 3.1415926535897932 + 0.5);

        my $mag_val = $star->{mag};
        if ($mag_val < -2.0) { $mag_val = -2.0; }
        if ($mag_val > 8.0) { $mag_val = 8.0; }
        my $mag_q = int(($mag_val + 2.0) * 255.0 / 10.0 + 0.5);

        print $ofh pack("V v v C C C C", $star->{hip}, $ra_q, $dec_q, $mag_q, 0, 0, 0);
    }
    
    # Write hash entries (uint32_t hash_key, uint16_t stars[4])
    foreach my $entry (@entry_records) {
        print $ofh pack("V v v v v", $entry->{hash_key}, $entry->{stars}->[0], $entry->{stars}->[1], $entry->{stars}->[2], $entry->{stars}->[3]);
    }
    
    close($ofh);
    print "Binary database compiled successfully!\n";
    print "  Total size: " . (-s $bin_db_path) / 1024 . " KB\n";
}

# Math helpers
sub asin {
    my $x = shift;
    return atan2($x, sqrt(1 - $x*$x));
}

sub acos_clamp {
    my $x = shift;
    $x = 1.0 if $x > 1.0;
    $x = -1.0 if $x < -1.0;
    return atan2(sqrt(1 - $x*$x), $x);
}

main();
