use strict;
use warnings;
use Test::More;
use File::Spec;
use JSON;

# Plan tests:
# Stage 1: Database verification (3 tests)
# Stage 2: Camera Simulation (3 tests)
# Stage 3: Centroid Extraction (3 tests)
# Stage 4: Star Identification (3 tests)
# Stage 5: Attitude Determination (4 tests)
# Stage 6: Celestial Position Fix (3 tests)
# Stage 7: EKF Trajectory Fusion (2 tests)
# Stage 8: Factor Graph Trajectory Smoothing (2 tests)
# Total tests: 23
plan tests => 23;

my $base_dir = File::Spec->catdir('data');
my $db_path = File::Spec->catfile($base_dir, 'starfix_db.bin');
my $mock_centroids_path = File::Spec->catfile($base_dir, 'mock_centroids.json');
my $mock_image_path = File::Spec->catfile($base_dir, 'mock_image.pgm');
my $detected_centroids_path = File::Spec->catfile($base_dir, 'detected_centroids.json');
my $identified_stars_path = File::Spec->catfile($base_dir, 'identified_stars.json');
my $attitude_path = File::Spec->catfile($base_dir, 'attitude.json');

# Helper to read JSON file
sub read_json {
    my $path = shift;
    open(my $fh, '<:encoding(UTF-8)', $path) or return undef;
    local $/;
    my $content = <$fh>;
    close($fh);
    return decode_json($content);
}

# --- Stage 1: Database verification ---
{
    ok(-f $db_path, "Database file exists");
    my $sz = -s $db_path;
    ok($sz > 10 * 1024 * 1024, "Database file size is larger than 10MB ($sz bytes)");
    
    open(my $fh, '<:raw', $db_path) or die "Cannot open database";
    my $buf;
    read($fh, $buf, 4);
    close($fh);
    
    ok($buf eq "SFIX", "Database has correct magic prefix SFIX");
}

# --- Stage 2: Camera Simulation ---
{
    my $output = `perl scripts/simulate_camera.pl 2>&1`;
    ok($? == 0, "Simulation step completed successfully");
    ok(-f $mock_centroids_path, "mock_centroids.json was generated");
    ok(-f $mock_image_path, "mock_image.pgm was generated");
}

# --- Stage 3: Centroid Extraction ---
{
    my $output = `./starfix_cli --centroid 2>&1`;
    ok($? == 0, "Centroiding step completed successfully");
    ok(-f $detected_centroids_path, "detected_centroids.json was generated");
    
    my $detected = read_json($detected_centroids_path);
    ok(defined $detected && @$detected > 0, "Centroid extraction detected stars");
}

# --- Stage 4: Star Identification ---
{
    my $output = `./starfix_cli --identify 2>&1`;
    ok($? == 0, "Star ID step completed successfully");
    ok(-f $identified_stars_path, "identified_stars.json was generated");
    
    my $identified = read_json($identified_stars_path);
    ok(defined $identified && @$identified >= 3, "TETRA Star ID matched at least 3 stars");
}

# --- Stage 5: Attitude Determination ---
{
    my $output = `./starfix_cli --estimate-pose 2>&1`;
    ok($? == 0, "Attitude estimation step completed successfully");
    ok(-f $attitude_path, "attitude.json was generated");
    
    my $att = read_json($attitude_path);
    ok(defined $att && exists $att->{ra_deg} && exists $att->{dec_deg} && exists $att->{roll_deg}, "Attitude contains pointing angles");
    my $err_ra = abs($att->{ra_deg} - 80.0);
    my $err_dec = abs($att->{dec_deg} - 37.7749);
    ok($err_ra < 0.1 && $err_dec < 0.1, "Pointing errors are within 0.1 degrees");
}

# --- Stage 6: Celestial Position Fix ---
{
    my $output = `./starfix_cli --solve-fix 2>&1`;
    ok($? == 0, "Position fix solver step completed successfully");
    
    ok($output =~ /Position Solver Error:\s+([\d\.]+)\s+Nautical Miles/, "Position solver error printed");
    if ($output =~ /Position Solver Error:\s+([\d\.]+)\s+Nautical Miles/) {
        ok($1 < 5.0, "Position fix solver error is under 5.0 nm ($1 nm)");
    } else {
        fail("Position fix accuracy check");
    }
}

# --- Stage 7: EKF Trajectory Fusion ---
{
    my $output = `./starfix_cli --fuse-ekf 2>&1`;
    ok($? == 0, "EKF trajectory fusion step completed successfully");
    ok($output =~ /Simulation Finished!/, "EKF completion message printed");
}

# --- Stage 8: Factor Graph Trajectory Smoothing ---
{
    my $output = `./starfix_cli --fuse-graph 2>&1`;
    ok($? == 0, "Factor Graph trajectory smoothing step completed successfully");
    ok($output =~ /Factor Graph Optimization Finished!/, "Factor Graph completion message printed");
}
