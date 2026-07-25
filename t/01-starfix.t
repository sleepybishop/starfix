use strict;
use warnings;
use Test::More;

# Plan tests:
# EKF (5 tests)
# Solver (3 tests)
# Centroiding (5 tests)
# Attitude (5 tests)
# Factor Graph (3 tests)
# Star ID (5 tests)
# Total tests: 26
plan tests => 26;

# Test 1: EKF implementation
{
    my $output = `./t/test_ekf 2>&1`;
    ok($? == 0, "EKF test runner executed successfully");
    
    ok($output =~ /Test 1: NULL pointer robustness.*Passed/, "EKF Test 1 (NULL pointers) passed");
    ok($output =~ /Test 2: Happy path propagation.*Passed/, "EKF Test 2 (DR propagation) passed");
    ok($output =~ /Test 3: Celestial correction.*Passed/, "EKF Test 3 (Celestial correction) passed");
    ok($output =~ /Test 4: Extreme noise correction limits.*Passed/, "EKF Test 4 (Covariance limits) passed");
}

# Test 2: Position Fix Solver (nanoqsp integration)
{
    my $output = `./t/test_solver 2>&1`;
    ok($? == 0, "Solver test runner executed successfully");
    
    ok($output =~ /Test 1: Argument validation.*Passed/, "Solver Test 1 (argument bounds) passed");
    ok($output =~ /Test 2: Happy path position fix solver.*Passed/, "Solver Test 2 (least squares convergence) passed");
}

# Test 3: Centroiding (Image processing)
{
    my $output = `./t/test_centroid 2>&1`;
    ok($? == 0, "Centroid test runner executed successfully");
    
    ok($output =~ /Test 1: Boundary check.*Passed/, "Centroid Test 1 (invalid bounds) passed");
    ok($output =~ /Test 2: Black image.*Passed/, "Centroid Test 2 (thresholding) passed");
    ok($output =~ /Test 3: Happy path multi-star detection.*Passed/, "Centroid Test 3 (weighted center of mass) passed");
    ok($output =~ /Test 4: Giant saturated blob rejection.*Passed/, "Centroid Test 4 (component size filter) passed");
}

# Test 4: Attitude Solver (Davenport's Q-method)
{
    my $output = `./t/test_attitude 2>&1`;
    ok($? == 0, "Attitude test runner executed successfully");
    
    ok($output =~ /Test 1: Boundary check.*Passed/, "Attitude Test 1 (boundary checks) passed");
    ok($output =~ /Test 2: Happy path attitude estimation.*Passed/, "Attitude Test 2 (Davenport eigenvalues) passed");
    ok($output =~ /Test 3: Collinear degenerate vectors.*Passed/, "Attitude Test 3 (degenerate inputs) passed");
    ok($output =~ /Test 6: Stellar Aberration Correction.*Passed/, "Attitude Test 4 (stellar aberration correction) passed");
}

# Test 5: Factor Graph Optimization (batch Gauss-Newton)
{
    my $output = `./t/test_graph 2>&1`;
    ok($? == 0, "Factor Graph test runner executed successfully");
    
    ok($output =~ /Test 1: Sizing and NULL validations.*Passed/, "Graph Test 1 (sizing bounds) passed");
    ok($output =~ /Test 2: Happy path factor graph optimization.*Passed/, "Graph Test 2 (Gauss-Newton optimization) passed");
}

# Test 6: Star Identification (TETRA Pattern Matching)
{
    my $output = `./t/test_identify 2>&1`;
    ok($? == 0, "Star ID test runner executed successfully");
    
    ok($output =~ /Test 1: Sizing and NULL validations.*Passed/s, "Star ID Test 1 (sizing bounds) passed");
    ok($output =~ /Test 2: Loading binary database.*Passed/s, "Star ID Test 2 (binary database loading) passed");
    ok($output =~ /Test 3: Lost-in-Space Star Identification.*Passed/s, "Star ID Test 3 (pattern search matching) passed");
    ok($output =~ /Test 3b: Tracking Mode with Attitude Hint.*Passed/s, "Star ID Test 4 (tracking mode with attitude hint) passed");
}
