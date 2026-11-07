#!/usr/bin/env perl
use strict;
use warnings;

my $width = 800;
my $height = 600;
my @img = (0) x ($width * $height);
my ($pole_x, $pole_y) = (400, 300);
my $pi = 3.14159265358979323846;

# draw concentric arcs
for my $radius (50, 100, 150, 200, 250, 300) {
    for (my $i = 0; $i < 100; $i++) {
        my $angle = ($pi / 2.0) * $i / 99.0;
        my $x = int($pole_x + $radius * cos($angle));
        my $y = int($pole_y + $radius * sin($angle));
        if ($x >= 0 && $x < $width && $y >= 0 && $y < $height) {
            $img[$y * $width + $x] = 255;
            # add some thickness
            $img[$y * $width + $x + 1] = 255 if $x + 1 < $width;
            $img[($y + 1) * $width + $x] = 255 if $y + 1 < $height;
        }
    }
}

open(my $fh, '>:raw', 'test_trails.pgm') or die "cannot write test_trails.pgm: $!";
print $fh "P5\n$width $height\n255\n", pack('C*', @img);
close($fh);
