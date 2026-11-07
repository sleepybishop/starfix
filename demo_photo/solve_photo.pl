#!/usr/bin/env perl
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use File::Spec;
use File::Temp qw(tempdir);

BEGIN {
    my $dir = dirname(abs_path(__FILE__));
    unshift @INC, $dir;
}
use photo_common;

my $PI = 3.14159265358979323846;

sub acos_clamp {
    my ($x) = @_;
    $x = 1.0 if $x > 1.0;
    $x = -1.0 if $x < -1.0;
    return atan2(sqrt(1.0 - $x * $x), $x);
}

sub tan_rad {
    my ($x) = @_;
    return sin($x) / cos($x);
}

sub get_camera_zenith {
    my ($tilts) = @_;
    $tilts ||= {};
    my $frame = $tilts->{camera_frame};
    if (defined $frame && $frame ne $photo_common::CAMERA_FRAME) {
        die "camera zenith must use +X right, +Y down, +Z forward\n";
    }
    my $zenith_camera = $tilts->{zenith_camera};
    if (defined $zenith_camera && ref($zenith_camera) eq 'ARRAY' && @$zenith_camera == 3) {
        for my $v (@$zenith_camera) {
            die "EXIF zenith_camera is not a unit vector\n"
                unless photo_common::is_finite_number($v);
        }
        my ($g_x, $g_y, $g_z) = map { $_ + 0.0 } @$zenith_camera;
        my $g_norm = sqrt($g_x * $g_x + $g_y * $g_y + $g_z * $g_z);
        if (!photo_common::is_finite_number($g_norm) || abs($g_norm - 1.0) > 1e-3) {
            die "EXIF zenith_camera is not a unit vector\n";
        }
        return [ $g_x / $g_norm, $g_y / $g_norm, $g_z / $g_norm ];
    }
    if (defined $zenith_camera) {
        die "zenith_camera must contain three components\n";
    }
    my $roll = defined $tilts->{roll_deg} ? $tilts->{roll_deg} : 0.0;
    my $pitch = defined $tilts->{pitch_deg} ? $tilts->{pitch_deg} : 0.0;
    if (!photo_common::is_finite_number($roll) || !photo_common::is_finite_number($pitch)) {
        die "non-finite legacy tilt component\n";
    }
    my $g_x = sin($roll * $PI / 180.0);
    my $g_y = sin($pitch * $PI / 180.0);
    my $g_z_sq = 1.0 - $g_x * $g_x - $g_y * $g_y;
    if ($g_z_sq < -1e-9) {
        die "legacy EXIF tilt components do not form a unit vector\n";
    }
    my $g_z = sqrt($g_z_sq > 0.0 ? $g_z_sq : 0.0);
    return [ $g_x, $g_y, $g_z ];
}

sub refract_camera_ray {
    my ($camera, $zenith) = @_;
    my $cos_theta = $camera->[0] * $zenith->[0]
        + $camera->[1] * $zenith->[1]
        + $camera->[2] * $zenith->[2];
    my $theta_t = acos_clamp($cos_theta);
    my $alt_t = $PI / 2.0 - $theta_t;
    return $camera if $theta_t <= 1e-6 || $alt_t <= 0.0;

    my $alt_deg = $alt_t * 180.0 / $PI;
    my $refraction;
    if ($alt_deg >= 15.0) {
        $refraction = 7.888888e-5 * 1010.0 / (283.0 * tan_rad($alt_t));
    } else {
        my $num = ((2e-5 * $alt_deg + 1.96e-2) * $alt_deg + 1.594e-1) * 1010.0;
        my $den = 283.0 * ((8.45e-2 * $alt_deg + 5.05e-1) * $alt_deg + 1.0);
        $refraction = ($num / $den) * $PI / 180.0;
    }
    my $theta_a = $theta_t - $refraction;
    $theta_a = 0.0 if $theta_a < 0.0;

    my @tangent = (
        $camera->[0] - $cos_theta * $zenith->[0],
        $camera->[1] - $cos_theta * $zenith->[1],
        $camera->[2] - $cos_theta * $zenith->[2],
    );
    my $norm = sqrt($tangent[0] * $tangent[0] + $tangent[1] * $tangent[1] + $tangent[2] * $tangent[2]);
    return $camera if $norm <= 1e-12;

    return [
        cos($theta_a) * $zenith->[0] + sin($theta_a) * ($tangent[0] / $norm),
        cos($theta_a) * $zenith->[1] + sin($theta_a) * ($tangent[1] / $norm),
        cos($theta_a) * $zenith->[2] + sin($theta_a) * ($tangent[2] / $norm),
    ];
}

sub solve_photo_overlay {
    my ($image_path, $zenith_camera) = @_;
    $image_path = abs_path($image_path) // die "cannot resolve $image_path\n";
    my $exif_dict = photo_common::load_exif($image_path);
    my $utc_dt = photo_common::parse_exif_utc($exif_dict);
    my (undef, $gha_aries) = photo_common::calculate_gha_aries($utc_dt);
    my $metadata = photo_common::get_user_comment($exif_dict->{'Exif'}) || {};

    if (defined $zenith_camera) {
        $metadata = { zenith_camera => [ @$zenith_camera ] };
    } elsif (!exists $metadata->{zenith_camera}
          && !exists $metadata->{roll_deg}
          && !exists $metadata->{pitch_deg}) {
        print "  No tilt measurement: assuming the image center points exactly at zenith\n";
    }
    my $zenith = get_camera_zenith($metadata);
    my $cli_bin = photo_common::get_cli_path();
    my $db_file = photo_common::get_db_path();
    my $iso_utc = photo_common::format_iso_utc($utc_dt);

    print "\nProcessing " . basename($image_path) . " at UTC $iso_utc\n";

    my $work_dir = tempdir('starfix-photo-XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $data_dir = File::Spec->catdir($work_dir, 'data');
    mkdir($data_dir) or die "cannot mkdir $data_dir: $!\n";
    symlink($db_file, File::Spec->catfile($data_dir, 'starfix_db.bin'))
        or die "cannot symlink starfix_db.bin: $!\n";

    my $plate = photo_common::solve_plate(
        $image_path, $exif_dict, $cli_bin, $work_dir, $utc_dt, $zenith
    );
    my $result = photo_common::run_command([
        $cli_bin, '--solve-photo',
        sprintf('%.9f', $gha_aries),
        sprintf('%.12f', $zenith->[0]),
        sprintf('%.12f', $zenith->[1]),
        sprintf('%.12f', $zenith->[2]),
        sprintf('%.9f', $plate->{fov_deg}),
        sprintf('%.9f', photo_common::decimal_year($utc_dt)),
        "$plate->{crop_size}",
    ], $work_dir);

    my $solution = photo_common::parse_solver_output($result);
    if (!photo_common::solution_is_acceptable($solution)) {
        die "photo position failed quality checks:\n$result->{stdout}$result->{stderr}";
    }
    my %position;
    for my $line (split(/\r?\n/, $result->{stdout})) {
        if ($line =~ /Solved Latitude\s*:\s*([+-]?\d+(?:\.\d+)?)/) {
            $position{latitude} = $1 + 0.0;
        } elsif ($line =~ /Solved Longitude\s*:\s*([+-]?\d+(?:\.\d+)?)/) {
            $position{longitude} = $1 + 0.0;
        }
    }
    if (scalar(keys %position) != 2) {
        die "solver did not return a complete position\n";
    }
    print $result->{stdout};

    my $width = $plate->{width};
    my $height = $plate->{height};
    my $size = $plate->{crop_size};
    my $left = int(($width - $size) / 2);
    my $top = int(($height - $size) / 2);
    my $focal = $size / (2.0 * tan_rad(($plate->{fov_deg} * $PI / 180.0) / 2.0));
    my $rotation = $solution->{rotation};

    my $im_cmd = photo_common::find_im_command();
    my @draw_cmd = (
        $im_cmd, $image_path, '-auto-orient',
        '-fill', 'none', '-stroke', 'lime', '-strokewidth', '3',
        '-draw', sprintf('rectangle %d,%d %d,%d', $left, $top, $left + $size, $top + $size),
    );

    for my $match (@{ $plate->{identified} }) {
        my $centroid = $plate->{centroids}[ $match->{centroid_idx} ];
        my $x = $left + $centroid->{u};
        my $y = $top + $centroid->{v};
        push @draw_cmd,
            '-fill', 'none', '-stroke', 'red', '-strokewidth', '2',
            '-draw', sprintf('circle %.2f,%.2f %.2f,%.2f', $x, $y, $x + 12.0, $y);

        my $ra = $match->{ra};
        my $dec = $match->{dec};
        my @celestial = (cos($dec) * cos($ra), cos($dec) * sin($ra), sin($dec));
        my @cam_raw = (
            $rotation->[0][0] * $celestial[0] + $rotation->[0][1] * $celestial[1] + $rotation->[0][2] * $celestial[2],
            $rotation->[1][0] * $celestial[0] + $rotation->[1][1] * $celestial[1] + $rotation->[1][2] * $celestial[2],
            $rotation->[2][0] * $celestial[0] + $rotation->[2][1] * $celestial[1] + $rotation->[2][2] * $celestial[2],
        );
        my $camera = refract_camera_ray(\@cam_raw, $zenith);
        next if $camera->[2] <= 0.0;

        my $px = $left + $size / 2.0 + $focal * ($camera->[0] / $camera->[2]);
        my $py = $top + $size / 2.0 + $focal * ($camera->[1] / $camera->[2]);
        push @draw_cmd,
            '-fill', 'cyan', '-stroke', 'none',
            '-draw', sprintf('circle %.2f,%.2f %.2f,%.2f', $px, $py, $px + 4.0, $py),
            '-fill', 'white', '-stroke', 'none', '-pointsize', '16',
            '-draw', sprintf("text %.2f,%.2f 'HIP %d'", $px + 16.0, $py + 5.0, $match->{hip});
    }

    my @text = (
        'STARFIX NAVIGATION',
        basename($image_path),
        "UTC: $iso_utc",
        sprintf('Camera zenith (+X right, +Y down, +Z forward): (%.6f, %.6f, %.6f)', @$zenith),
        sprintf('Latitude: %.6f  Longitude: %.6f', $position{latitude}, $position{longitude}),
        sprintf('Crop: %dpx  RMS: %.3fpx  Inliers: %d/%d',
                $size, $solution->{rms}, $solution->{inliers}, $solution->{matches}),
    );
    my $box_right = ($width - 20 < 1050) ? ($width - 20) : 1050;
    push @draw_cmd,
        '-fill', 'rgba(0,0,0,0.82)', '-stroke', 'none',
        '-draw', sprintf('rectangle 20,20 %d,180', $box_right),
        '-fill', 'white', '-stroke', 'none', '-pointsize', '16';
    for (my $i = 0; $i < @text; $i++) {
        (my $line = $text[$i]) =~ s/'/\\'/g;
        push @draw_cmd, '-draw', sprintf("text 30,%d '%s'", 44 + 23 * $i, $line);
    }

    (my $out_jpeg = $image_path) =~ s/\.[^.\/\\]+\z//;
    $out_jpeg .= '_overlay.jpg';
    push @draw_cmd, '-quality', '85', $out_jpeg;
    system(@draw_cmd) == 0 or die "failed to render overlay image $out_jpeg\n";
    print "Saved $out_jpeg\n";
    return $out_jpeg;
}

sub main {
    my (@argv) = @_;
    my @images;
    my $zenith;

    for (my $i = 0; $i < @argv; $i++) {
        if ($argv[$i] eq '--zenith') {
            if ($i + 3 >= @argv) {
                warn "error: --zenith requires 3 numeric arguments\n";
                return 2;
            }
            my @z = @argv[$i + 1 .. $i + 3];
            for my $v (@z) {
                if (!photo_common::is_finite_number($v)) {
                    warn "error: invalid --zenith value: $v\n";
                    return 2;
                }
            }
            $zenith = [ map { $_ + 0.0 } @z ];
            $i += 3;
        } elsif ($argv[$i] =~ /\A-/) {
            warn "error: unknown option $argv[$i]\n";
            return 2;
        } else {
            push @images, $argv[$i];
        }
    }

    if (!@images) {
        my $search_dir = -d 'demo_photo' ? 'demo_photo' : '.';
        @images = grep { !/_overlay\.jpg\z/ } sort glob("$search_dir/*.jpg");
    }
    if (!@images) {
        warn "error: no JPEG images found\n";
        return 2;
    }

    my $failed = 0;
    for my $path (@images) {
        eval { solve_photo_overlay($path, $zenith); 1 } or do {
            my $err = $@ || 'unknown error';
            chomp $err;
            warn "Failed $path: $err\n";
            $failed++;
        };
    }
    return $failed ? 2 : 0;
}

if (!caller()) {
    exit(main(@ARGV));
}

1;
