use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempdir);
use JSON qw(decode_json encode_json);
use POSIX qw(floor);
use Test::More;

my $root_dir = dirname(dirname(abs_path(__FILE__)));
BEGIN {
    my $root = dirname(dirname(abs_path(__FILE__)));
    unshift @INC, File::Spec->catdir($root, 'demo_photo');
}
use photo_common;
require File::Spec->catfile($root_dir, 'demo_photo', 'solve_photo.pl');

my $PI = 3.14159265358979323846;

# test 1: FOV is invariant under portrait/landscape rotation
{
    my $exif = { 'Exif' => { $photo_common::TAG_FOCAL_LENGTH_35MM => 24 } };
    my $portrait = photo_common::estimate_crop_fov($exif, 3000, 4000);
    my $landscape = photo_common::estimate_crop_fov($exif, 4000, 3000);
    ok(abs($portrait - $landscape) < 1e-12, 'FOV is invariant under portrait/landscape rotation');
}

# test 2: real photo geometry and mirror rejection
{
    my $fixture_path = File::Spec->catfile($root_dir, 't', 'fixtures', 'photo_centroids.json');
    open(my $fh, '<:encoding(UTF-8)', $fixture_path) or die "cannot read $fixture_path: $!";
    local $/;
    my $fixture = decode_json(<$fh>);
    close($fh);

    my ($latitude, $longitude) = @{ $fixture->{known_position} };
    for my $case (@{ $fixture->{photos} }) {
        my $directory = tempdir('starfix-test-photo-XXXXXX', TMPDIR => 1, CLEANUP => 1);
        my $data_dir = File::Spec->catdir($directory, 'data');
        mkdir($data_dir) or die "mkdir $data_dir: $!";
        symlink(
            File::Spec->catfile($root_dir, 'data', 'starfix_db.bin'),
            File::Spec->catfile($data_dir, 'starfix_db.bin')
        ) or die "symlink starfix_db.bin: $!";

        my $centroids = $case->{centroids};
        my $centroids_file = File::Spec->catfile($data_dir, 'detected_centroids.json');
        open(my $cfh, '>:encoding(UTF-8)', $centroids_file) or die $!;
        print $cfh encode_json($centroids);
        close($cfh);

        my $utc = photo_common::parse_iso_utc($case->{utc});
        my (undef, $gha) = photo_common::calculate_gha_aries($utc);
        my $year = photo_common::decimal_year($utc);
        my $cli = File::Spec->catfile($root_dir, 'starfix_cli');

        my $res = photo_common::run_command(
            [ $cli, '--identify', "$case->{size}", "$case->{size}", "$case->{fov_deg}", "$year" ],
            $directory
        );
        is($res->{returncode}, 0, "$case->{name}: --identify succeeded") or diag($res->{stdout});

        my @fit_args = ($cli, '--fit-photo', '0', '0', '1', "$case->{fov_deg}", "$year", "$case->{size}");
        $res = photo_common::run_command(\@fit_args, $directory);
        is($res->{returncode}, 0, "$case->{name}: --fit-photo succeeded") or diag($res->{stdout});

        my $solution = photo_common::parse_solver_output($res);
        ok(photo_common::solution_is_acceptable($solution), "$case->{name}: solution is acceptable")
            or diag($res->{stdout});
        cmp_ok($solution->{inliers}, '>=', 5, "$case->{name}: at least 5 inliers");

        my $rotation = $solution->{rotation};
        my ($ref_lat, $ref_lon) = @{ $case->{reference_center_lat_lon} };
        my $dec = $ref_lat * $PI / 180.0;
        my $ra = ($ref_lon + $gha) * $PI / 180.0;
        my @expected = (cos($dec) * cos($ra), cos($dec) * sin($ra), sin($dec));
        my $dot = $rotation->[2][0] * $expected[0]
            + $rotation->[2][1] * $expected[1]
            + $rotation->[2][2] * $expected[2];
        my $error = acos_clamp($dot) * 180.0 / $PI;
        cmp_ok($error, '<', 0.08, "$case->{name}: center pointing error < 0.08 deg");

        my $z_dec = $latitude * $PI / 180.0;
        my $z_ra = ($longitude + $gha) * $PI / 180.0;
        my @zenith = (cos($z_dec) * cos($z_ra), cos($z_dec) * sin($z_ra), sin($z_dec));
        my $z_camera = $rotation->[2][0] * $zenith[0]
            + $rotation->[2][1] * $zenith[1]
            + $rotation->[2][2] * $zenith[2];
        my $tilt = acos_clamp($z_camera) * 180.0 / $PI;
        cmp_ok(abs($tilt - $case->{reference_tilt_deg}), '<=', 0.08, "$case->{name}: tilt matches reference");

        # the same matched stars reflected vertically cannot be explained by a rotation
        my @mirrored = map {
            +{ %$_, v => $case->{size} - $_->{v} }
        } @$centroids;
        open(my $mfh, '>:encoding(UTF-8)', $centroids_file) or die $!;
        print $mfh encode_json(\@mirrored);
        close($mfh);

        $res = photo_common::run_command(\@fit_args, $directory);
        isnt($res->{returncode}, 0, "$case->{name}: reflected stars rejected");
        unlike($res->{stdout}, qr/ROTATION_MATRIX:/, "$case->{name}: reflected stars emit no rotation matrix");
    }
}

# test 3: solver output requires complete quality data
{
    my $result = {
        stdout     => "RMS_ERROR: 0.25 pixels\nINLIERS: 4/5\nROTATION_MATRIX: 1 0 0 0 1 0 0 0 1\n",
        returncode => 0,
    };
    my $solution = photo_common::parse_solver_output($result);
    ok(photo_common::solution_is_acceptable($solution), 'complete quality output accepted');
    $solution->{inliers} = 2;
    ok(!photo_common::solution_is_acceptable($solution), 'too few inliers rejected');
    $solution->{inliers} = 4;
    $solution->{rms} = 9**9**9;
    ok(!photo_common::solution_is_acceptable($solution), 'infinite RMS rejected');
}

# test 4: EXIF time uses offset and subseconds
{
    my $exif = {
        'Exif' => {
            $photo_common::TAG_DATETIME_ORIGINAL    => '2024:07:09 22:30:15',
            $photo_common::TAG_SUBSEC_TIME_ORIGINAL => '25',
            $photo_common::TAG_OFFSET_TIME_ORIGINAL => '-07:00',
        },
        '0th' => {},
    };
    my $utc = photo_common::parse_exif_utc($exif);
    is(photo_common::format_iso_utc($utc), '2024-07-10T05:30:15.250000', 'EXIF time applies offset and subseconds');
}

# test 5: negative zenith hemisphere round-trips
{
    my @expected = (0.2, -0.3, -sqrt(0.87));
    my $actual = get_camera_zenith({ zenith_camera => \@expected });
    for my $i (0 .. 2) {
        ok(abs($expected[$i] - $actual->[$i]) < 1e-12, "zenith component $i round-trips");
    }
}

# test 6: EXIF orientation is applied before crop
{
    my $directory = tempdir('starfix-orient-XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $source_path = File::Spec->catfile($directory, 'portrait.jpg');
    my $output_path = File::Spec->catfile($directory, 'crop.pgm');
    my $im_cmd = photo_common::find_im_command();
    system($im_cmd, '-size', '600x800', 'xc:black', $source_path) == 0
        or die "failed to create test JPEG\n";
    my $exif = {
        '0th'  => { $photo_common::TAG_ORIENTATION => 6 },
        'Exif' => {},
    };
    photo_common::write_jpeg_exif($source_path, $exif);
    my ($width, $height) = photo_common::crop_center_to_pgm($source_path, $output_path);
    is($width, 800, 'oriented width is 800');
    is($height, 600, 'oriented height is 600');

    open(my $pfh, '<:raw', $output_path) or die "cannot open $output_path: $!";
    my $magic = <$pfh>;
    my $dims = <$pfh>;
    close($pfh);
    chomp $dims if defined $dims;
    is($dims, "$photo_common::CROP_SIZE $photo_common::CROP_SIZE", 'cropped PGM has CROP_SIZE x CROP_SIZE');
}

# test 7: sparse center uses wider crop
{
    my $fixture_path = File::Spec->catfile($root_dir, 't', 'fixtures', 'photo_centroids.json');
    open(my $fh, '<:encoding(UTF-8)', $fixture_path) or die "cannot read $fixture_path: $!";
    local $/;
    my $fixture = decode_json(<$fh>);
    close($fh);
    my $case = $fixture->{photos}[0];

    my $directory = tempdir('starfix-sparse-XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $size = $case->{size};
    my @pixels = (0) x ($size * $size);
    for my $c (@{ $case->{centroids} }) {
        my $cx = int($c->{u} + 0.5);
        my $cy = int($c->{v} + 0.5);
        for my $dy (-2 .. 2) {
            for my $dx (-2 .. 2) {
                next if $dx * $dx + $dy * $dy > 4;
                my $px = $cx + $dx;
                my $py = $cy + $dy;
                if ($px >= 0 && $px < $size && $py >= 0 && $py < $size) {
                    $pixels[$py * $size + $px] = 220;
                }
            }
        }
    }

    my $pgm_path = File::Spec->catfile($directory, 'sparse_src.pgm');
    open(my $pfh, '>:raw', $pgm_path) or die "cannot write $pgm_path: $!";
    print $pfh "P5\n$size $size\n255\n", pack('C*', @pixels);
    close($pfh);

    my $jpg_path = File::Spec->catfile($directory, 'sparse.jpg');
    my $im_cmd = photo_common::find_im_command();
    system($im_cmd, $pgm_path, '-quality', '100', $jpg_path) == 0
        or die "failed to convert sparse PGM to JPEG\n";

    my $utc = photo_common::parse_iso_utc($case->{utc});
    my $fov_rad = $case->{fov_deg} * $PI / 180.0;
    my $tan_half = sin($fov_rad / 2.0) / cos($fov_rad / 2.0);
    my $focal_mm = int(
        ($size / (2.0 * $tan_half)) * sqrt(36.0 * 36.0 + 24.0 * 24.0) / sqrt(2.0 * $size * $size) + 0.5
    );
    my $date_str = sprintf(
        '%04d:%02d:%02d %02d:%02d:%02d',
        $utc->{year}, $utc->{month}, $utc->{day},
        $utc->{hour}, $utc->{minute}, $utc->{second}
    );
    my $exif = {
        '0th'  => {},
        'Exif' => {
            $photo_common::TAG_DATETIME_ORIGINAL    => $date_str,
            $photo_common::TAG_OFFSET_TIME_ORIGINAL => '+00:00',
            $photo_common::TAG_FOCAL_LENGTH_35MM    => $focal_mm,
        },
    };
    photo_common::write_jpeg_exif($jpg_path, $exif);

    open(my $orig_fh, '<:raw', $jpg_path) or die $!;
    my $original_bytes = do { local $/; <$orig_fh> };
    close($orig_fh);

    my $data_dir = File::Spec->catdir($directory, 'data');
    mkdir($data_dir) or die "mkdir $data_dir: $!";
    symlink(
        File::Spec->catfile($root_dir, 'data', 'starfix_db.bin'),
        File::Spec->catfile($data_dir, 'starfix_db.bin')
    ) or die "symlink starfix_db.bin: $!";

    my $cli = File::Spec->catfile($root_dir, 'starfix_cli');
    my $plate = photo_common::solve_plate($jpg_path, $exif, $cli, $directory, $utc);
    is($plate->{crop_size}, 840, 'sparse center selects 840px crop');

    open(my $after_fh, '<:raw', $jpg_path) or die $!;
    my $after_bytes = do { local $/; <$after_fh> };
    close($after_fh);
    is($after_bytes, $original_bytes, 'solve_plate does not modify source JPEG');

    my $solve_script = File::Spec->catfile($root_dir, 'demo_photo', 'solve_photo.pl');
    my $res = photo_common::run_command(
        [ $^X, $solve_script, '--zenith', '0', '0', '1', $jpg_path ],
        $root_dir
    );
    is($res->{returncode}, 0, 'solve_photo.pl succeeded on sparse JPEG') or diag($res->{stdout});

    my ($solved_lat, $solved_lon);
    for my $line (split(/\r?\n/, $res->{stdout})) {
        $solved_lat = $1 + 0.0 if $line =~ /Solved Latitude\s*:\s*([+-]?\d+(?:\.\d+)?)/;
        $solved_lon = $1 + 0.0 if $line =~ /Solved Longitude\s*:\s*([+-]?\d+(?:\.\d+)?)/;
    }
    my ($expected_lat, $expected_lon) = @{ $case->{reference_center_lat_lon} };
    cmp_ok(abs($solved_lat - $expected_lat), '<=', 0.08, 'sparse JPEG solved latitude within 0.08 deg');
    cmp_ok(abs($solved_lon - $expected_lon), '<=', 0.08, 'sparse JPEG solved longitude within 0.08 deg');

    my $overlay_path = File::Spec->catfile($directory, 'sparse_overlay.jpg');
    ok(-f $overlay_path, 'sparse_overlay.jpg was created');

    open(my $final_fh, '<:raw', $jpg_path) or die $!;
    my $final_bytes = do { local $/; <$final_fh> };
    close($final_fh);
    is($final_bytes, $original_bytes, 'solve_photo.pl does not modify source JPEG');
}

done_testing();
