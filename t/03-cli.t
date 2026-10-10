use strict;
use warnings;
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
done_testing();
