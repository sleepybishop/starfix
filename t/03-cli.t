use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Copy qw(copy);
use File::Temp qw(tempdir);
use Test::More;

my $dir = tempdir(CLEANUP => 1);
is(system('./t/test_cli_io', "$dir/records.json"), 0, 'Flat record readers, bounds, and matrix layout');
for my $command ('--fit-photo', '--solve-photo') {
    my $output = `./starfix_cli $command 0 0 2>&1`;
    ok($? != 0 && $output =~ /Usage:/, "$command rejects missing arguments");
    $output = `./starfix_cli $command 0 0 1 12 2026 560 1 1 2>&1`;
    ok($? != 0 && $output =~ /Usage:/, "$command rejects extra arguments");
}
my $output = `./starfix_cli --identify 1024 1024 12 nan 2>&1`;
ok($? != 0, 'Identification rejects a nonfinite observation year');

my $root = getcwd();
my $cli = abs_path('starfix_cli');
sub run_cli {
    my $pid = open(my $fh, '-|');
    die "fork: $!" unless defined $pid;
    if (!$pid) {
        open(STDERR, '>&', STDOUT) or die "redirect stderr: $!";
        exec {$cli} $cli, @_;
        die "exec: $!";
    }
    local $/;
    my $text = <$fh>;
    close($fh);
    return ($?, $text);
}

for my $args (
    ['--centroid', 'garbage'], ['--centroid', '256'], ['--centroid', '-1'],
    ['--centroid', '1.5'], ['--centroid', '35junk'], ['--centroid', ''],
    ['--centroid', '999999999999999999999999'], ['--centroid', '35', 'extra'],
    ['--identify', '0'], ['--identify', '1024.5'], ['--identify', '1e100'],
    ['--identify', '1024', '1024', '12junk'],
    ['--identify', '1024', '1024', '12', '2026', 'extra'],
    ['--fit-photo', 'garbage', '0', '1'], ['--fit-photo', '0', '0', '1e9999'],
    ['--fit-photo', '0', '0', '1', '1e-9999'],
    ['--solve-photo', 'garbage', '0', '0', '1'],
) {
    my ($status, $text) = run_cli(@$args);
    ok($status != 0 && $text =~ /Invalid numeric argument|Usage:/,
       "Reject invalid arguments: @$args") or diag($text);
}

my $work = tempdir(TMPDIR => 1, CLEANUP => 1);
mkdir "$work/data" or die $!;
open(my $pgm, '>:raw', "$work/data/mock_image.pgm") or die $!;
print {$pgm} "P5\n3 3\n255\n", pack('C*', 0, 0, 0, 0, 100, 0, 0, 0, 0);
close($pgm) or die $!;
copy('data/starfix_db.bin', "$work/data/starfix_db.bin") or die $!;
chdir($work) or die $!;
for my $threshold (0, 255) {
    my ($status, $text) = run_cli('--centroid', $threshold);
    is($status, 0, "Accept threshold boundary $threshold") or diag($text);
}
unlink 'data/detected_centroids.json' or die $!;
mkdir 'data/detected_centroids.json' or die $!;
my ($status, $text) = run_cli('--centroid');
ok($status != 0 && $text =~ /Could not write/, 'Centroid CLI reports output open failure');
rmdir 'data/detected_centroids.json' or die $!;
copy("$root/data/detected_centroids.json", 'data/detected_centroids.json') or die $!;
mkdir 'data/identified_stars.json' or die $!;
($status, $text) = run_cli('--identify');
ok($status != 0 && $text =~ /Could not write/ && $text !~ /Saved \d/,
   'Identify CLI reports output open failure without claiming it saved results') or diag($text);
rmdir 'data/identified_stars.json' or die $!;
SKIP: {
    skip '/dev/full is unavailable', 2 unless -c '/dev/full';
    symlink '/dev/full', 'data/identified_stars.json' or die $!;
    ($status, $text) = run_cli('--identify');
    ok($status != 0 && $text =~ /Could not write/, 'Identify CLI reports buffered write failure');
    unlink 'data/detected_centroids.json' or die $!;
    symlink '/dev/full', 'data/detected_centroids.json' or die $!;
    ($status, $text) = run_cli('--centroid');
    ok($status != 0 && $text =~ /Could not write/, 'Centroid CLI reports buffered write failure');
}
unlink 'data/detected_centroids.json' or die $!;
unlink 'data/identified_stars.json' if -l 'data/identified_stars.json';
copy("$root/data/detected_centroids.json", 'data/detected_centroids.json') or die $!;
copy("$root/data/identified_stars.json", 'data/identified_stars.json') or die $!;
mkdir 'child' or die $!;
mkdir 'child/data' or die $!;
chdir 'child' or die $!;
open(my $bad, '>', 'data/detected_centroids.json') or die $!;
print {$bad} '[{"u":"garbage","v":0}]';
close($bad) or die $!;
($status, $text) = run_cli('--fit-photo', 0, 0, 1);
ok($status != 0 && $text =~ /Could not load centroids/,
   'Malformed local centroids do not fall back to parent data');
unlink 'data/detected_centroids.json' or die $!;
open($bad, '>', 'data/identified_stars.json') or die $!;
print {$bad} '[{"centroid_idx":"garbage","catalog_idx":0}]';
close($bad) or die $!;
($status, $text) = run_cli('--fit-photo', 0, 0, 1);
ok($status != 0 && $text =~ /found: -4/,
   'Malformed local matches do not fall back to parent data');
chdir($root) or die $!;
done_testing();
