#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use JSON;

sub get_camera_rotation_matrix {
    my ($ra_rad, $dec_rad, $roll_rad) = @_;
    
    # z axis points toward target ra/dec
    my $z_cam = [
        cos($dec_rad) * cos($ra_rad),
        cos($dec_rad) * sin($ra_rad),
        sin($dec_rad)
    ];
    
    # x axis points east (zero roll horizontal axis)
    my $x_zero = [
        -sin($ra_rad),
        cos($ra_rad),
        0.0
    ];
    
    # y axis points north (zero roll vertical axis) via cross product
    my $y_zero = [
        $z_cam->[1]*$x_zero->[2] - $z_cam->[2]*$x_zero->[1],
        $z_cam->[2]*$x_zero->[0] - $z_cam->[0]*$x_zero->[2],
        $z_cam->[0]*$x_zero->[1] - $z_cam->[1]*$x_zero->[0]
    ];
    
    # apply roll angle around z axis
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

# box-muller transform to generate standard normal variables
sub rand_normal {
    my ($mean, $std) = @_;
    my $u1 = rand();
    my $u2 = rand();
    while ($u1 <= 1e-15) { $u1 = rand(); } # avoid log(0)
    my $z = sqrt(-2.0 * log($u1)) * cos(2.0 * 3.1415926535897932 * $u2);
    return $mean + $z * $std;
}

sub main {
    my $script_dir = File::Spec->rel2abs(File::Spec->path());
    # find project base dir
    my $base_dir = File::Spec->catdir(File::Spec->splitdir(File::Spec->rel2abs($0)));
    my @dirs = File::Spec->splitdir($base_dir);
    pop @dirs; # pop filename
    pop @dirs; # pop scripts dir
    $base_dir = File::Spec->catdir(@dirs);
    
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
    
    my @stars;
    my $catalog_buf;
    my $catalog_bytes = $num_stars * 12;
    read($fh, $catalog_buf, $catalog_bytes) == $catalog_bytes or die "Short read on star catalog";
    
    my $PI = 3.1415926535897932;
    for (my $i = 0; $i < $num_stars; $i++) {
        my $star_entry = substr($catalog_buf, $i * 12, 12);
        my ($hip, $ra_q, $dec_q, $mag_q, $p1, $p2, $p3, $p4) = unpack("V v v C C C C", $star_entry);
        
        my $ra = $ra_q * (2.0 * $PI) / 65535.0;
        my $dec = $dec_q * $PI / 65535.0 - $PI / 2.0;
        my $mag = $mag_q * 10.0 / 255.0 - 2.0;
        
        push @stars, {
            id   => $i,
            hip  => $hip,
            name => "Unknown",
            ra   => $ra,
            dec  => $dec,
            mag  => $mag
        };
    }
    close($fh);
    
    print "Loaded " . scalar(@stars) . " catalog stars.\n";
    
    my $ra_deg = 80.0;
    my $dec_deg = 37.7749;
    my $roll_deg = 0.0;
    my $fov_deg = 12.0;
    my $width = 1024;
    my $height = 1024;
    
    my $ra_rad = $ra_deg * 3.1415926535897932 / 180.0;
    my $dec_rad = $dec_deg * 3.1415926535897932 / 180.0;
    my $roll_rad = $roll_deg * 3.1415926535897932 / 180.0;
    my $fov_rad = $fov_deg * 3.1415926535897932 / 180.0;
    
    print "Simulating camera pointing:\n";
    print "  RA = $ra_deg deg, Dec = $dec_deg deg, Roll = $roll_deg deg\n";
    print "  FOV = $fov_deg deg, Sensor = ${width}x${height} pixels\n";
    
    my ($x_cam, $y_cam, $z_cam) = get_camera_rotation_matrix($ra_rad, $dec_rad, $roll_rad);
    my $focal_length = $width / (2.0 * tax($fov_rad / 2.0));
    
    # helper for tan since Perl has no core tan
    sub tax {
        my $val = shift;
        return sin($val) / cos($val);
    }
    
    my $year = 1900.0 + (gmtime)[5] + (gmtime)[7] / 365.25;
    my $T = ($year - 2000.0) / 100.0;
    my $PI = 3.1415926535897932;
    my $zeta = 2306.2181 * $T * $PI / (3600.0 * 180.0);
    my $z = 2306.2181 * $T * $PI / (3600.0 * 180.0);
    my $theta = 2004.3109 * $T * $PI / (3600.0 * 180.0);
    
    my $c_zeta = cos($zeta); my $s_zeta = sin($zeta);
    my $c_z = cos($z); my $s_z = sin($z);
    my $c_th = cos($theta); my $s_th = sin($theta);
    
    my $P00 = $c_zeta * $c_th * $c_z - $s_zeta * $s_z;
    my $P01 = -$s_zeta * $c_th * $c_z - $c_zeta * $s_z;
    my $P02 = -$s_th * $c_z;
    my $P10 = $c_zeta * $c_th * $s_z + $s_zeta * $c_z;
    my $P11 = -$s_zeta * $c_th * $s_z + $c_zeta * $c_z;
    my $P12 = -$s_th * $s_z;
    my $P20 = $c_zeta * $s_th;
    my $P21 = -$s_zeta * $s_th;
    my $P22 = $c_th;
    
    my @visible_stars;
    my $idx = 0;
    foreach my $star (@stars) {
        my $ra_s = $star->{ra};
        my $dec_s = $star->{dec};
        
        my $v0_0 = cos($dec_s) * cos($ra_s);
        my $v0_1 = cos($dec_s) * sin($ra_s);
        my $v0_2 = sin($dec_s);
        
        my $v_star = [
            $P00 * $v0_0 + $P01 * $v0_1 + $P02 * $v0_2,
            $P10 * $v0_0 + $P11 * $v0_1 + $P12 * $v0_2,
            $P20 * $v0_0 + $P21 * $v0_1 + $P22 * $v0_2
        ];
        
        # dot product for rotation projection
        my $x_c = $x_cam->[0]*$v_star->[0] + $x_cam->[1]*$v_star->[1] + $x_cam->[2]*$v_star->[2];
        my $y_c = $y_cam->[0]*$v_star->[0] + $y_cam->[1]*$v_star->[1] + $y_cam->[2]*$v_star->[2];
        my $z_c = $z_cam->[0]*$v_star->[0] + $z_cam->[1]*$v_star->[1] + $z_cam->[2]*$v_star->[2];
        
        if ($z_c > 0) {
            my $x_sensor = $x_c / $z_c;
            my $y_sensor = $y_c / $z_c;
            
            my $u = $width / 2.0 + $focal_length * $x_sensor;
            my $v = $height / 2.0 - $focal_length * $y_sensor;
            
            if ($u >= 0 && $u < $width && $v >= 0 && $v < $height) {
                push @visible_stars, {
                    db_index => $idx,
                    id => $star->{id},
                    hip => $star->{hip},
                    name => $star->{name} || "Unknown",
                    mag => $star->{mag},
                    u_true => $u,
                    v_true => $v
                };
            }
        }
        $idx++;
    }
    
    print "Found " . scalar(@visible_stars) . " visible stars in camera view.\n";
    
    # save mock centroids to JSON
    open(my $cfh, '>', $centroids_path) or die "Cannot write centroids: $!";
    print $cfh encode_json(\@visible_stars);
    close($cfh);
    print "Saved true visible star list to $centroids_path\n";
    
    # render mock image
    print "Rendering mock star field image...\n";
    
    # allocate flat 1D image array
    my @img = (0) x ($width * $height);
    
    foreach my $star (@visible_stars) {
        my $u = $star->{u_true};
        my $v = $star->{v_true};
        my $mag = $star->{mag};
        
        my $r = 6.0 - $mag;
        $r = int($r + 0.5);
        if ($r < 1) { $r = 1; }
        
        my $peak_intensity = int(255.0 * (10.0 ** (-0.4 * ($mag - 1.0))) + 0.5);
        if ($peak_intensity < 80) { $peak_intensity = 80; }
        if ($peak_intensity > 255) { $peak_intensity = 255; }
        
        for (my $y_offset = -$r; $y_offset <= $r; $y_offset++) {
            for (my $x_offset = -$r; $x_offset <= $r; $x_offset++) {
                my $dist_sq = $x_offset**2 + $y_offset**2;
                if ($dist_sq <= $r**2) {
                    my $pu = int($u + $x_offset + 0.5);
                    my $pv = int($v + $y_offset + 0.5);
                    if ($pu >= 0 && $pu < $width && $pv >= 0 && $pv < $height) {
                        my $sigma = ($r > 1) ? ($r / 2.0) : 0.5;
                        my $factor = exp(-$dist_sq / (2.0 * $sigma**2));
                        my $val = int($peak_intensity * $factor + 0.5);
                        
                        my $px_idx = $pv * $width + $pu;
                        $img[$px_idx] += $val;
                        if ($img[$px_idx] > 255) { $img[$px_idx] = 255; }
                    }
                }
            }
        }
    }
    
    # add background and Gaussian noise (reduce std from 4 to 0.5 for sub-pixel precision)
    srand(42); # fix seed for reproducibility
    for (my $i = 0; $i < $width * $height; $i++) {
        my $noise = rand_normal(15, 0.5);
        my $val = int($img[$i] + $noise + 0.5);
        if ($val < 0) { $val = 0; }
        if ($val > 255) { $val = 255; }
        $img[$i] = $val;
    }
    
    # save as PGM (P5 binary format)
    open(my $ifh, '>:raw', $image_path) or die "Cannot write image: $!";
    print $ifh "P5\n$width $height\n255\n";
    my $packed = pack("C*", @img);
    print $ifh $packed;
    close($ifh);
    print "Saved rendered noisy star field image to $image_path\n";
}

main();
