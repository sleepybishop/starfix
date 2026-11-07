use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename qw(dirname);
use File::Copy qw(copy);
use File::Find qw(find);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

my $root_dir = dirname(dirname(abs_path(__FILE__)));

sub write_text {
    my ($path, $content) = @_;
    open(my $fh, '>', $path) or die "cannot write $path: $!";
    print $fh $content;
    close($fh);
}

sub read_text {
    my ($path) = @_;
    open(my $fh, '<', $path) or die "cannot read $path: $!";
    local $/;
    my $text = <$fh>;
    close($fh);
    return $text;
}

sub read_lines {
    my ($path) = @_;
    return [] unless -e $path;
    my $text = read_text($path);
    my @lines = split(/\r?\n/, $text);
    return \@lines;
}

sub set_all_mtime {
    my ($dir, $ts) = @_;
    find(sub {
        if (-f $_) {
            utime($ts, $ts, $_) or die "utime $_: $!";
        }
    }, $dir);
}

{
    my $temp = tempdir('starfix-build-XXXXXX', TMPDIR => 1, CLEANUP => 1);
    my $prototype = File::Spec->catdir($temp, 'prototype');
    mkdir($prototype) or die "mkdir $prototype: $!";
    copy(
        File::Spec->catfile($root_dir, 'prototype', 'Makefile'),
        File::Spec->catfile($prototype, 'Makefile')
    ) or die "copy prototype/Makefile: $!";

    write_text(
        File::Spec->catfile($temp, 'Makefile'),
        ".PHONY: record\n"
        . "record:\n\t\@echo checked >> checks.log\n"
        . "libstarfix.a: core.c | record\n\tcp \$< \$\@\n"
        . "libstarfix_navigation.a: navigation.c | record\n\tcp \$< \$\@\n"
    );

    for my $name ('core.c', 'navigation.c') {
        write_text(File::Spec->catfile($temp, $name), $name);
    }
    for my $name ('main', 'camera', 'bno055', 'sim', 'camera_mock', 'bno055_mock', 'trail_preprocessor') {
        write_text(File::Spec->catfile($prototype, "$name.c"), $name);
    }

    my $compiler = File::Spec->catfile($temp, 'compiler.pl');
    write_text(
        $compiler,
        "use strict;\nuse warnings;\nuse File::Basename qw(basename);\n"
        . "my \$out;\n"
        . "for (my \$i = 0; \$i < \@ARGV; \$i++) {\n"
        . "    if (\$ARGV[\$i] eq '-o') { \$out = \$ARGV[\$i + 1]; last; }\n"
        . "}\n"
        . "open(my \$fh, '>', \$out) or die \$!;\n"
        . "print \$fh 'built';\n"
        . "close(\$fh);\n"
        . "open(my \$log, '>>', 'compiler.log') or die \$!;\n"
        . "print \$log basename(\$out) . \"\\n\";\n"
        . "close(\$log);\n"
    );

    my $build = sub {
        local %ENV = %ENV;
        delete $ENV{MAKEFLAGS};
        delete $ENV{MFLAGS};
        my $cc = "$^X $compiler";
        my $output = `make -j8 -C "$prototype" "CC=$cc" 2>&1`;
        is($?, 0, 'standalone prototype make succeeded') or diag($output);
    };

    my $checks = File::Spec->catfile($temp, 'checks.log');
    my $compilations = File::Spec->catfile($prototype, 'compiler.log');

    $build->();
    is_deeply(read_lines($checks), ['checked'], 'initial build checks root libraries once');
    my $initial_compilations = read_lines($compilations);
    ok((grep { $_ eq 'starfix_prototype' } @$initial_compilations) > 0, 'built starfix_prototype');
    ok((grep { $_ eq 'starfix_prototype_sim' } @$initial_compilations) > 0, 'built starfix_prototype_sim');

    # fixed old timestamps make update checks independent of clock resolution
    my $baseline = time() - 3600;
    set_all_mtime($temp, $baseline);
    my $before = read_lines($compilations);
    $build->();
    is_deeply(read_lines($checks), ['checked', 'checked'], 'second build checks root once');
    is_deeply(read_lines($compilations), $before, 'unchanged build does not relink');

    my $index = 3;
    for my $pair (['core.c', 'libstarfix.a'], ['navigation.c', 'libstarfix_navigation.a']) {
        my ($source, $archive) = @$pair;
        set_all_mtime($temp, $baseline);
        my $changed = File::Spec->catfile($temp, $source);
        write_text($changed, "updated $source");
        utime($baseline + 1, $baseline + 1, $changed) or die "utime $changed: $!";

        $before = read_lines($compilations);
        $build->();
        is(
            read_text(File::Spec->catfile($temp, $archive)),
            read_text($changed),
            "$archive refreshed from $source"
        );
        is_deeply(
            read_lines($checks),
            [ ('checked') x $index ],
            "root check count after updating $source"
        );
        my $after = read_lines($compilations);
        my @relinked = sort @{$after}[ scalar(@$before) .. $#$after ];
        is_deeply(
            \@relinked,
            ['starfix_prototype', 'starfix_prototype_sim'],
            "relinked both prototype targets when $archive changed"
        );
        $index++;
    }

    $before = read_lines($compilations);
    $build->();
    is_deeply(read_lines($checks), [ ('checked') x 5 ], 'final unchanged build checks root once');
    is_deeply(read_lines($compilations), $before, 'final unchanged build does not relink');
}

done_testing();
