#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use JSON;

my $PI = 3.14159265358979323846;

sub get_camera_rotation_matrix {
    my ($ra_rad, $dec_rad, $roll_rad) = @_;

    my $z_cam = [
        cos($dec_rad) * cos($ra_rad),
        cos($dec_rad) * sin($ra_rad),
        sin($dec_rad)
    ];
    my $x_zero = [ -sin($ra_rad), cos($ra_rad), 0.0 ];
    my $y_zero = [
        $z_cam->[1] * $x_zero->[2] - $z_cam->[2] * $x_zero->[1],
        $z_cam->[2] * $x_zero->[0] - $z_cam->[0] * $x_zero->[2],
        $z_cam->[0] * $x_zero->[1] - $z_cam->[1] * $x_zero->[0]
    ];
    my $x_cam = [
        cos($roll_rad) * $x_zero->[0] + sin($roll_rad) * $y_zero->[0],
        cos($roll_rad) * $x_zero->[1] + sin($roll_rad) * $y_zero->[1],
        cos($roll_rad) * $x_zero->[2] + sin($roll_rad) * $y_zero->[2]
    ];
    my $y_cam = [
        -sin($roll_rad) * $x_zero->[0] + cos($roll_rad) * $y_zero->[0],
        -sin($roll_rad) * $x_zero->[1] + cos($roll_rad) * $y_zero->[1],
        -sin($roll_rad) * $x_zero->[2] + cos($roll_rad) * $y_zero->[2]
    ];
    return ($x_cam, $y_cam, $z_cam);
}

sub rand_normal {
    my ($mean, $std) = @_;
    my $u1 = rand();
    my $u2 = rand();
    $u1 = rand() while $u1 <= 1e-15;
    my $z = sqrt(-2.0 * log($u1)) * cos(2.0 * $PI * $u2);
    return $mean + $z * $std;
}

sub tan_rad {
    my ($val) = @_;
    return sin($val) / cos($val);
}

sub astrometry_frame {
    my ($year) = @_;
    my $mjd = 51544.5 + ($year - 2000.0) * 365.25;
    my $T = ($mjd - 51544.5) / 36525.0;
    my $arcsec_rad = $PI / (180.0 * 3600.0);

    my $zeta = $arcsec_rad * (2.650545 + $T * (2306.083227 + $T * (0.2988499 + $T * 0.01801828)));
    my $z = $arcsec_rad * (-2.650545 + $T * (2306.077181 + $T * (1.0927348 + $T * 0.01826837)));
    my $theta = $arcsec_rad * ($T * (2004.191903 - $T * (0.4294934 + $T * 0.04182264)));

    my ($cz, $sz) = (cos($zeta), sin($zeta));
    my ($cza, $sza) = (cos($z), sin($z));
    my ($ct, $st) = (cos($theta), sin($theta));
    my @prec = (
        [ $cz * $ct * $cza - $sz * $sza, -$sz * $ct * $cza - $cz * $sza, -$st * $cza ],
        [ $cz * $ct * $sza + $sz * $cza, -$sz * $ct * $sza + $cz * $cza, -$st * $sza ],
        [ $cz * $st,                     -$sz * $st,                     $ct ],
    );

    my $eps = (23.4392911 + $T * (-46.8150 / 3600.0)) * $PI / 180.0;
    my $omega = (125.04452 - 1934.136261 * $T) * $PI / 180.0;
    my $Lsun = (280.4665 + 36000.7698 * $T) * $PI / 180.0;
    my $Lmoon = (218.3165 + 481267.8813 * $T) * $PI / 180.0;
    my $dpsi = (-17.20 * sin($omega) - 1.32 * sin(2.0 * $Lsun)
                - 0.23 * sin(2.0 * $Lmoon) + 0.21 * sin(2.0 * $omega)) * $arcsec_rad;
    my $deps = (9.20 * cos($omega) + 0.57 * cos(2.0 * $Lsun)
                + 0.10 * cos(2.0 * $Lmoon) - 0.09 * cos(2.0 * $omega)) * $arcsec_rad;

    my ($se, $ce) = (sin($eps), cos($eps));
    my ($sp, $cp) = (sin($dpsi), cos($dpsi));
    my ($sede, $cede) = (sin($eps + $deps), cos($eps + $deps));
    my @nut = (
        [ $cp,          -$sp * $ce,                 -$sp * $se ],
        [ $cede * $sp,  $cede * $cp * $ce + $sede * $se, $cede * $cp * $se - $sede * $ce ],
        [ $sede * $sp,  $sede * $cp * $ce - $cede * $se, $sede * $cp * $se + $cede * $ce ],
    );

    my @P;
    for my $i (0 .. 2) {
        for my $j (0 .. 2) {
            $P[$i][$j] = $nut[$i][0] * $prec[0][$j]
                       + $nut[$i][1] * $prec[1][$j]
                       + $nut[$i][2] * $prec[2][$j];
        }
    }

    my $L = 2.0 * $PI * (0.27908 + 100.00214 * $T);
    my @aberr = (-0.994e-4 * sin($L), 0.912e-4 * cos($L), 0.395e-4 * cos($L));
    return (\@P, \@aberr);
}

sub main {
    my @dirs = File::Spec->splitdir(File::Spec->rel2abs($0));
    pop @dirs;
    pop @dirs;
    my $base_dir = File::Spec->catdir(@dirs);

    my $db_path = File::Spec->catfile($base_dir, "data", "starfix_db.bin");
    my $centroids_path = File::Spec->catfile($base_dir, "data", "mock_centroids.json");
    my $image_path = File::Spec->catfile($base_dir, "data", "mock_image.pgm");

    print "Loading database from $db_path...\n";
    open(my $fh, '<:raw', $db_path) or die "Cannot open database: $!";

    my $magic;
    read($fh, $magic, 4) == 4 or die "Short read on magic prefix";
    die "Invalid database magic prefix: $magic" unless $magic eq "SFIX";

    my $header_buf;
    read($fh, $header_buf, 12) == 12 or die "Short read on header";
    my ($num_stars, $num_entries, $bin_factor) = unpack("V V V", $header_buf);

    my $catalog_bytes = $num_stars * 18;
    my $catalog_buf;
    read($fh, $catalog_buf, $catalog_bytes) == $catalog_bytes or die "Short read on star catalog";
    close($fh);

    my $mas_to_rad = $PI / (180.0 * 3600000.0);
    my @stars;
    for (my $i = 0; $i < $num_stars; $i++) {
        my $entry = substr($catalog_buf, $i * 18, 18);
        my ($hip, $ra_q, $dec_q, $mag_q, $pm_ra_q, $pm_dec_q) = unpack("l< V V s< s< s<", $entry);
        push @stars, {
            id     => $i,
            hip    => $hip,
            name   => "Unknown",
            ra     => $ra_q * (2.0 * $PI) / 4294967295.0,
            dec    => $dec_q * $PI / 4294967295.0 - $PI / 2.0,
            mag    => $mag_q * 1e-3,
            pm_ra  => $pm_ra_q * 0.2 * $mas_to_rad,
            pm_dec => $pm_dec_q * 0.2 * $mas_to_rad,
        };
    }

    print "Loaded " . scalar(@stars) . " catalog stars.\n";

    my $ra_deg = 80.0;
    my $dec_deg = 37.7749;
    my $roll_deg = 0.0;
    my $fov_deg = 12.0;
    my $width = 1024;
    my $height = 1024;

    my $ra_rad = $ra_deg * $PI / 180.0;
    my $dec_rad = $dec_deg * $PI / 180.0;
    my $roll_rad = $roll_deg * $PI / 180.0;
    my $fov_rad = $fov_deg * $PI / 180.0;

    print "Simulating camera pointing:\n";
    print "  RA = $ra_deg deg, Dec = $dec_deg deg, Roll = $roll_deg deg\n";
    print "  FOV = $fov_deg deg, Sensor = ${width}x${height} pixels\n";

    my ($x_cam, $y_cam, $z_cam) = get_camera_rotation_matrix($ra_rad, $dec_rad, $roll_rad);
    my $focal_length = $width / (2.0 * tan_rad($fov_rad / 2.0));

    my $year = 1900.0 + (gmtime)[5] + (gmtime)[7] / 365.25;
    my $dt_yr = $year - 2000.0;
    my ($P, $aberr) = astrometry_frame($year);

    my @visible_stars;
    my $idx = 0;
    foreach my $star (@stars) {
        my ($ra0, $dec0) = ($star->{ra}, $star->{dec});
        my ($ca, $sa) = (cos($ra0), sin($ra0));
        my ($cd, $sd) = (cos($dec0), sin($dec0));
        my $pm_a = $star->{pm_ra} * $dt_yr;
        my $pm_d = $star->{pm_dec} * $dt_yr;
        my @r0 = (
            $cd * $ca - $sa * $pm_a - $sd * $ca * $pm_d,
            $cd * $sa + $ca * $pm_a - $sd * $sa * $pm_d,
            $sd + $cd * $pm_d,
        );
        my $r0_norm = sqrt($r0[0] * $r0[0] + $r0[1] * $r0[1] + $r0[2] * $r0[2]);
        @r0 = map { $_ / $r0_norm } @r0 if $r0_norm > 0.0;

        my @v_star = map {
            $P->[$_][0] * $r0[0] + $P->[$_][1] * $r0[1] + $P->[$_][2] * $r0[2] + $aberr->[$_]
        } 0 .. 2;
        my $v_norm = sqrt($v_star[0] * $v_star[0] + $v_star[1] * $v_star[1] + $v_star[2] * $v_star[2]);
        @v_star = map { $_ / $v_norm } @v_star if $v_norm > 0.0;

        my $x_c = $x_cam->[0] * $v_star[0] + $x_cam->[1] * $v_star[1] + $x_cam->[2] * $v_star[2];
        my $y_c = $y_cam->[0] * $v_star[0] + $y_cam->[1] * $v_star[1] + $y_cam->[2] * $v_star[2];
        my $z_c = $z_cam->[0] * $v_star[0] + $z_cam->[1] * $v_star[1] + $z_cam->[2] * $v_star[2];

        if ($z_c > 0.0) {
            my $u = $width / 2.0 + $focal_length * ($x_c / $z_c);
            my $v = $height / 2.0 + $focal_length * ($y_c / $z_c);
            if ($u >= 0.0 && $u < $width && $v >= 0.0 && $v < $height) {
                push @visible_stars, {
                    db_index => $idx,
                    id       => $star->{id},
                    hip      => $star->{hip},
                    name     => $star->{name},
                    mag      => $star->{mag},
                    u_true   => $u,
                    v_true   => $v,
                };
            }
        }
        $idx++;
    }

    @visible_stars = sort { $a->{mag} <=> $b->{mag} } @visible_stars;
    print "Found " . scalar(@visible_stars) . " visible stars in camera view.\n";

    open(my $cfh, '>', $centroids_path) or die "Cannot write centroids: $!";
    print $cfh encode_json(\@visible_stars);
    close($cfh);
    print "Saved true visible star list to $centroids_path\n";

    print "Rendering mock star field image...\n";
    my @img = (0) x ($width * $height);

    foreach my $star (@visible_stars) {
        my ($u, $v, $mag) = ($star->{u_true}, $star->{v_true}, $star->{mag});
        my $r = 3;
        my $sigma = 1.0;

        my $peak_intensity = 240.0 * (10.0 ** (-0.20 * ($mag - 1.0)));
        $peak_intensity = 45.0 if $peak_intensity < 45.0;
        $peak_intensity = 240.0 if $peak_intensity > 240.0;

        my $cu = int($u + 0.5);
        my $cv = int($v + 0.5);
        for (my $y_offset = -$r; $y_offset <= $r; $y_offset++) {
            for (my $x_offset = -$r; $x_offset <= $r; $x_offset++) {
                my $pu = $cu + $x_offset;
                my $pv = $cv + $y_offset;
                if ($pu >= 0 && $pu < $width && $pv >= 0 && $pv < $height) {
                    my $dist_sq = ($pu - $u)**2 + ($pv - $v)**2;
                    my $val = int($peak_intensity * exp(-$dist_sq / (2.0 * $sigma**2)) + 0.5);
                    my $px_idx = $pv * $width + $pu;
                    $img[$px_idx] += $val;
                    $img[$px_idx] = 255 if $img[$px_idx] > 255;
                }
            }
        }
    }

    srand(42);
    for (my $i = 0; $i < $width * $height; $i++) {
        my $val = int($img[$i] + rand_normal(15, 0.5) + 0.5);
        $val = 0 if $val < 0;
        $val = 255 if $val > 255;
        $img[$i] = $val;
    }

    open(my $ifh, '>:raw', $image_path) or die "Cannot write image: $!";
    print $ifh "P5\n$width $height\n255\n";
    print $ifh pack("C*", @img);
    close($ifh);
    print "Saved rendered noisy star field image to $image_path\n";
}

main();
