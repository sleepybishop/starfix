package photo_common;
use strict;
use warnings;

use Cwd qw(abs_path);
use File::Spec;
use JSON qw(decode_json encode_json);
use POSIX qw(floor);
use Time::Local qw(timegm);

our $CROP_SIZE = 560;
our $MAX_RMS_PIXELS = 3.0;
our $MIN_INLIERS = 3;
our $CAMERA_FRAME = 'x_right_y_down_z_forward';

# standard TIFF/EXIF tag constants
our $TAG_ORIENTATION = 274;
our $TAG_DATETIME = 306;
our $TAG_EXIF_IFD = 34665;
our $TAG_DATETIME_ORIGINAL = 36867;
our $TAG_OFFSET_TIME = 36880;
our $TAG_OFFSET_TIME_ORIGINAL = 36881;
our $TAG_OFFSET_TIME_DIGITIZED = 36882;
our $TAG_USER_COMMENT = 37510;
our $TAG_SUBSEC_TIME_ORIGINAL = 37521;
our $TAG_FOCAL_LENGTH_35MM = 41989;

my $PI = 3.14159265358979323846;

sub is_finite_number {
    my ($val) = @_;
    return 0 unless defined $val && !ref($val);
    return 0 unless $val =~ /\A[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?\z/;
    my $num = $val + 0.0;
    return ($num == $num && $num + 0.0 != 9**9**9 && $num + 0.0 != -(9**9**9)) ? 1 : 0;
}

sub find_im_command {
    return $ENV{IMTOOL} if $ENV{IMTOOL};
    for my $tool ('magick', 'convert') {
        for my $dir (File::Spec->path()) {
            my $candidate = File::Spec->catfile($dir, $tool);
            return $candidate if -f $candidate && -x $candidate;
        }
    }
    die "ImageMagick (magick or convert) was not found in PATH\n";
}

sub make_utc {
    my ($year, $month, $day, $hour, $minute, $second, $microsecond) = @_;
    $microsecond ||= 0;
    my $epoch = timegm($second, $minute, $hour, $day, $month - 1, $year);
    return {
        year        => $year + 0,
        month       => $month + 0,
        day         => $day + 0,
        hour        => $hour + 0,
        minute      => $minute + 0,
        second      => $second + 0,
        microsecond => $microsecond + 0,
        epoch       => $epoch + $microsecond / 1_000_000.0,
    };
}

sub epoch_us_to_utc {
    my ($epoch_us) = @_;
    my $whole_sec = floor($epoch_us / 1_000_000);
    my $us = int(sprintf('%.0f', $epoch_us - $whole_sec * 1_000_000));
    if ($us >= 1_000_000) {
        $whole_sec += 1;
        $us -= 1_000_000;
    }
    my ($sec, $min, $hour, $mday, $mon, $year) = gmtime($whole_sec);
    return make_utc($year + 1900, $mon + 1, $mday, $hour, $min, $sec, $us);
}

sub parse_iso_utc {
    my ($text) = @_;
    die "invalid ISO timestamp\n" unless defined $text;
    if ($text =~ /\A(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2})(?:\.(\d+))?Z?\z/) {
        my ($y, $m, $d, $hh, $mm, $ss, $frac) = ($1, $2, $3, $4, $5, $6, $7);
        my $us = 0;
        if (defined $frac && length($frac) > 0) {
            my $digits = substr($frac . '000000', 0, 6);
            $us = int($digits);
        }
        return make_utc($y, $m, $d, $hh, $mm, $ss, $us);
    }
    die "invalid ISO timestamp: $text\n";
}

sub format_iso_utc {
    my ($utc) = @_;
    my $base = sprintf(
        '%04d-%02d-%02dT%02d:%02d:%02d',
        $utc->{year}, $utc->{month}, $utc->{day},
        $utc->{hour}, $utc->{minute}, $utc->{second}
    );
    if (($utc->{microsecond} || 0) > 0) {
        return sprintf('%s.%06d', $base, $utc->{microsecond});
    }
    return $base;
}

sub calculate_gha_aries {
    my ($utc_time) = @_;
    my $year = $utc_time->{year};
    my $month = $utc_time->{month};
    my $day = $utc_time->{day};
    my $hour = $utc_time->{hour};
    my $minute = $utc_time->{minute};
    my $second = $utc_time->{second} + ($utc_time->{microsecond} || 0) / 1_000_000.0;

    my $a = floor((14 - $month) / 12);
    my $y = $year + 4800 - $a;
    my $m = $month + 12 * $a - 3;
    my $jdn = $day + floor((153 * $m + 2) / 5) + 365 * $y
        + floor($y / 4) - floor($y / 100) + floor($y / 400) - 32045;
    my $jd = $jdn + ($hour - 12) / 24.0 + $minute / 1440.0 + $second / 86400.0;

    my $t = ($jd - 2451545.0) / 36525.0;
    my $gmst = 280.46061837
        + 360.98564736629 * ($jd - 2451545.0)
        + 0.000387933 * $t * $t
        - ($t * $t * $t) / 38710000.0;

    # equation of the equinoxes (GAST - GMST) matching true-of-date catalog frame
    my $eps = (23.4392911 + $t * (-46.8150 / 3600.0)) * $PI / 180.0;
    my $omega = (125.04452 - 1934.136261 * $t) * $PI / 180.0;
    my $lsun = (280.4665 + 36000.7698 * $t) * $PI / 180.0;
    my $lmoon = (218.3165 + 481267.8813 * $t) * $PI / 180.0;
    my $dpsi_deg = (-17.20 * sin($omega) - 1.32 * sin(2.0 * $lsun)
                    - 0.23 * sin(2.0 * $lmoon) + 0.21 * sin(2.0 * $omega)) / 3600.0;
    my $deps = (9.20 * cos($omega) + 0.57 * cos(2.0 * $lsun)
                + 0.10 * cos(2.0 * $lmoon) - 0.09 * cos(2.0 * $omega)) * $PI / (180.0 * 3600.0);
    my $gast = $gmst + $dpsi_deg * cos($eps + $deps);
    $gast = $gast - 360.0 * floor($gast / 360.0);
    return ($jd, $gast);
}

sub decimal_year {
    my ($utc_time) = @_;
    my $year = $utc_time->{year};
    my $start = timegm(0, 0, 0, 1, 0, $year);
    my $end = timegm(0, 0, 0, 1, 0, $year + 1);
    my $current = defined $utc_time->{epoch}
        ? $utc_time->{epoch}
        : timegm(
            $utc_time->{second}, $utc_time->{minute}, $utc_time->{hour},
            $utc_time->{day}, $utc_time->{month} - 1, $year
        ) + ($utc_time->{microsecond} || 0) / 1_000_000.0;
    return $year + ($current - $start) / ($end - $start);
}

sub load_exif {
    my ($image_path) = @_;
    open(my $fh, '<:raw', $image_path) or die "cannot open $image_path: $!\n";
    local $/;
    my $data = <$fh>;
    close($fh);

    die "not a JPEG file: $image_path\n"
        unless defined $data && length($data) >= 4 && substr($data, 0, 2) eq "\xFF\xD8";

    my %result = ('0th' => {}, 'Exif' => {}, 'GPS' => {}, '1st' => {});
    my $pos = 2;
    while ($pos + 4 <= length($data)) {
        last if substr($data, $pos, 1) ne "\xFF";
        $pos++ while $pos < length($data) && substr($data, $pos, 1) eq "\xFF";
        last if $pos >= length($data);
        my $marker = ord(substr($data, $pos++, 1));
        last if $marker == 0xDA || $marker == 0xD9;
        last if $pos + 2 > length($data);
        my $seg_len = unpack('n', substr($data, $pos, 2));
        last if $seg_len < 2 || $pos + $seg_len > length($data);
        my $seg = substr($data, $pos + 2, $seg_len - 2);

        if ($marker == 0xE1 && length($seg) >= 14 && substr($seg, 0, 6) eq "Exif\x00\x00") {
            my $tiff = substr($seg, 6);
            my $byte_order = substr($tiff, 0, 2);
            if ($byte_order eq 'II' || $byte_order eq 'MM') {
                my ($u16, $u32) = $byte_order eq 'II' ? ('v', 'V') : ('n', 'N');
                my $ifd0_off = unpack($u32, substr($tiff, 4, 4));
                my $parse_ifd = sub {
                    my ($off) = @_;
                    my %tags;
                    return \%tags if $off + 2 > length($tiff);
                    my $count = unpack($u16, substr($tiff, $off, 2));
                    my %type_sizes = (1 => 1, 2 => 1, 3 => 2, 4 => 4, 5 => 8, 7 => 1, 9 => 4, 10 => 8);
                    for (my $i = 0; $i < $count; $i++) {
                        my $entry_off = $off + 2 + $i * 12;
                        last if $entry_off + 12 > length($tiff);
                        my $entry = substr($tiff, $entry_off, 12);
                        my ($tag, $type, $num) = unpack("$u16$u16$u32", substr($entry, 0, 8));
                        my $unit = $type_sizes{$type} || 1;
                        my $total = $unit * $num;
                        my $raw;
                        if ($total <= 4) {
                            $raw = substr($entry, 8, $total);
                        } else {
                            my $val_off = unpack($u32, substr($entry, 8, 4));
                            next if $val_off + $total > length($tiff);
                            $raw = substr($tiff, $val_off, $total);
                        }
                        if ($type == 2) {
                            $raw =~ s/\x00+\z//;
                            $tags{$tag} = $raw;
                        } elsif ($type == 3 && $num == 1) {
                            $tags{$tag} = unpack($u16, $raw);
                        } elsif ($type == 4 && $num == 1) {
                            $tags{$tag} = unpack($u32, $raw);
                        } else {
                            $tags{$tag} = $raw;
                        }
                    }
                    return \%tags;
                };
                $result{'0th'} = $parse_ifd->($ifd0_off);
                if (my $exif_off = $result{'0th'}{$TAG_EXIF_IFD}) {
                    $result{'Exif'} = $parse_ifd->($exif_off);
                }
            }
        } elsif (($marker >= 0xC0 && $marker <= 0xC3)
              || ($marker >= 0xC5 && $marker <= 0xC7)
              || ($marker >= 0xC9 && $marker <= 0xCB)
              || ($marker >= 0xCD && $marker <= 0xCF)) {
            if (length($seg) >= 5) {
                my ($h, $w) = unpack('nn', substr($seg, 1, 4));
                $result{_raw_width} = $w;
                $result{_raw_height} = $h;
            }
        }
        $pos += $seg_len;
    }
    return \%result;
}

sub pack_ifd_section {
    my ($tags_ref, $base_offset) = @_;
    my @keys = sort { $a <=> $b } keys %$tags_ref;
    my $count = scalar(@keys);
    my $entries_bytes = '';
    my $extra_bytes = '';
    my $data_offset = $base_offset + 2 + $count * 12 + 4;

    for my $tag (@keys) {
        my $val = $tags_ref->{$tag};
        my ($type, $num, $raw);
        if ($tag == $TAG_USER_COMMENT) {
            $type = 7;
            $raw = $val =~ /\AASCII\x00\x00\x00/ ? $val : ("ASCII\x00\x00\x00" . $val);
            $num = length($raw);
        } elsif ($tag == $TAG_EXIF_IFD) {
            $type = 4;
            $num = 1;
            $raw = pack('V', $val);
        } elsif (defined $val && $val =~ /\A\d+\z/ && $val <= 65535) {
            $type = 3;
            $num = 1;
            $raw = pack('v', $val);
        } else {
            $type = 2;
            $raw = $val . "\x00";
            $num = length($raw);
        }
        my $val_field;
        if (length($raw) <= 4) {
            $val_field = substr($raw . "\x00\x00\x00\x00", 0, 4);
        } else {
            $val_field = pack('V', $data_offset + length($extra_bytes));
            $extra_bytes .= $raw;
            $extra_bytes .= "\x00" if length($extra_bytes) % 2;
        }
        $entries_bytes .= pack('vvV', $tag, $type, $num) . $val_field;
    }
    return pack('v', $count) . $entries_bytes . pack('V', 0) . $extra_bytes;
}

sub strip_jpeg_exif_segments {
    my ($data) = @_;
    my $pos = 2;
    while ($pos + 4 <= length($data)) {
        last if substr($data, $pos, 1) ne "\xFF";
        my $marker_start = $pos;
        $pos++ while $pos < length($data) && substr($data, $pos, 1) eq "\xFF";
        last if $pos >= length($data);
        my $marker = ord(substr($data, $pos++, 1));
        last if $marker == 0xDA || $marker == 0xD9;
        if ($marker == 0x01 || ($marker >= 0xD0 && $marker <= 0xD7)) {
            next;
        }
        last if $pos + 2 > length($data);
        my $seg_len = unpack('n', substr($data, $pos, 2));
        last if $seg_len < 2 || $pos + $seg_len > length($data);
        if ($marker == 0xE1 && $seg_len >= 8 && substr($data, $pos + 2, 6) eq "Exif\x00\x00") {
            substr($data, $marker_start, ($pos + $seg_len) - $marker_start, '');
            $pos = $marker_start;
        } else {
            $pos += $seg_len;
        }
    }
    return $data;
}

sub write_jpeg_exif {
    my ($image_path, $exif_dict) = @_;
    open(my $fh, '<:raw', $image_path) or die "cannot open $image_path: $!\n";
    local $/;
    my $data = <$fh>;
    close($fh);
    die "not a JPEG file: $image_path\n" unless substr($data, 0, 2) eq "\xFF\xD8";
    $data = strip_jpeg_exif_segments($data);

    my %ifd0 = %{ $exif_dict->{'0th'} || {} };
    my %exif = %{ $exif_dict->{'Exif'} || {} };

    # calculate IFD0 size to determine where Exif IFD starts
    if (keys %exif) {
        $ifd0{$TAG_EXIF_IFD} = 0;
        my $dummy_ifd0 = pack_ifd_section(\%ifd0, 8);
        my $exif_offset = 8 + length($dummy_ifd0);
        $ifd0{$TAG_EXIF_IFD} = $exif_offset;
        my $ifd0_bytes = pack_ifd_section(\%ifd0, 8);
        my $exif_bytes = pack_ifd_section(\%exif, $exif_offset);
        my $tiff = 'II' . pack('vV', 42, 8) . $ifd0_bytes . $exif_bytes;
        my $app1 = "\xFF\xE1" . pack('n', length($tiff) + 8) . "Exif\x00\x00" . $tiff;
        substr($data, 2, 0, $app1);
    } else {
        my $ifd0_bytes = pack_ifd_section(\%ifd0, 8);
        my $tiff = 'II' . pack('vV', 42, 8) . $ifd0_bytes;
        my $app1 = "\xFF\xE1" . pack('n', length($tiff) + 8) . "Exif\x00\x00" . $tiff;
        substr($data, 2, 0, $app1);
    }

    open(my $out, '>:raw', $image_path) or die "cannot write $image_path: $!\n";
    print $out $data;
    close($out);
}

sub get_oriented_dimensions {
    my ($image_path, $exif_dict) = @_;
    $exif_dict ||= load_exif($image_path);
    my $w = $exif_dict->{_raw_width};
    my $h = $exif_dict->{_raw_height};
    if (!defined $w || !defined $h) {
        my $full = load_exif($image_path);
        $w = $full->{_raw_width};
        $h = $full->{_raw_height};
        $exif_dict->{'0th'}{$TAG_ORIENTATION} //= $full->{'0th'}{$TAG_ORIENTATION};
    }
    die "could not determine JPEG dimensions for $image_path\n" unless defined $w && defined $h;
    my $orient = $exif_dict->{'0th'}{$TAG_ORIENTATION} || 1;
    if ($orient >= 5 && $orient <= 8) {
        ($w, $h) = ($h, $w);
    }
    return ($w, $h);
}

sub get_user_comment {
    my ($exif_ifd) = @_;
    return undef unless $exif_ifd && ref($exif_ifd) eq 'HASH';
    my $raw = $exif_ifd->{$TAG_USER_COMMENT};
    return undef unless defined $raw;
    my $text = $raw;
    if (length($text) >= 8 && substr($text, 0, 5) eq 'ASCII') {
        $text = substr($text, 8);
    }
    $text =~ s/\x00+\z//;
    my $decoded = eval { decode_json($text) };
    return (defined $decoded && ref($decoded) eq 'HASH') ? $decoded : undef;
}

sub parse_exif_utc {
    my ($exif_dict) = @_;
    my $exif_ifd = $exif_dict->{'Exif'} || {};
    my $ifd0 = $exif_dict->{'0th'} || {};

    my $date_value = $exif_ifd->{$TAG_DATETIME_ORIGINAL};
    my @offset_tags = ($TAG_OFFSET_TIME_ORIGINAL, $TAG_OFFSET_TIME, $TAG_OFFSET_TIME_DIGITIZED);
    if (!defined $date_value) {
        $date_value = $ifd0->{$TAG_DATETIME};
        @offset_tags = ($TAG_OFFSET_TIME, $TAG_OFFSET_TIME_ORIGINAL, $TAG_OFFSET_TIME_DIGITIZED);
    }
    die "no EXIF timestamp\n" unless defined $date_value;

    die "invalid EXIF timestamp: $date_value\n"
        unless $date_value =~ /\A(\d{4}):(\d{2}):(\d{2}) (\d{2}):(\d{2}):(\d{2})\z/;
    my ($year, $mon, $day, $hour, $min, $sec) = ($1, $2, $3, $4, $5, $6);

    my $microsecond = 0;
    my $subsec = $exif_ifd->{$TAG_SUBSEC_TIME_ORIGINAL};
    if (defined $subsec && $subsec ne '') {
        (my $digits = $subsec) =~ s/\D+//g;
        $digits = substr($digits, 0, 6);
        if (length($digits) > 0) {
            $microsecond = int(substr($digits . '000000', 0, 6));
        }
    }

    my $offset_hours;
    for my $tag (@offset_tags) {
        my $value = $exif_ifd->{$tag};
        next unless defined $value && length($value) == 6;
        next unless $value =~ /\A([+-])(\d{2}):(\d{2})\z/;
        my ($sign_ch, $hh, $mm) = ($1, int($2), int($3));
        next if $hh > 14 || $mm > 59 || ($hh == 14 && $mm != 0);
        my $sign = $sign_ch eq '+' ? 1.0 : -1.0;
        $offset_hours = $sign * ($hh + $mm / 60.0);
        last;
    }

    if (!defined $offset_hours) {
        my $configured = $ENV{STARFIX_UTC_OFFSET_HOURS};
        if (defined $configured && is_finite_number($configured)) {
            $offset_hours = $configured + 0.0;
        }
    }
    if (!defined $offset_hours || !is_finite_number($offset_hours)
        || $offset_hours < -14.0 || $offset_hours > 14.0) {
        die "EXIF timestamp has no valid UTC offset; set STARFIX_UTC_OFFSET_HOURS explicitly\n";
    }

    my $local_epoch = timegm($sec, $min, $hour, $day, $mon - 1, $year);
    my $epoch_us = $local_epoch * 1_000_000 + $microsecond
        - int(sprintf('%.0f', $offset_hours * 3_600_000_000.0));
    return epoch_us_to_utc($epoch_us);
}

sub crop_center_to_pgm {
    my ($image_path, $output_path, $crop_size) = @_;
    $crop_size //= $CROP_SIZE;
    my ($width, $height) = get_oriented_dimensions($image_path);
    if ($crop_size < 1 || $width < $crop_size || $height < $crop_size) {
        die "image is smaller than the required ${crop_size}x${crop_size} crop\n";
    }
    my $left = int(($width - $crop_size) / 2);
    my $top = int(($height - $crop_size) / 2);
    my $im_cmd = find_im_command();
    my @cmd = (
        $im_cmd, $image_path,
        '-auto-orient',
        '-crop', sprintf('%dx%d+%d+%d', $crop_size, $crop_size, $left, $top),
        '+repage',
        '-grayscale', 'Rec601Luma',
        '-depth', '8',
        "pgm:$output_path",
    );
    system(@cmd) == 0 or die "failed to crop $image_path to PGM\n";
    return ($width, $height);
}

sub get_cli_path {
    for my $candidate ('./starfix_cli', '../starfix_cli') {
        if (-f $candidate && -x $candidate) {
            return abs_path($candidate);
        }
    }
    for my $dir (File::Spec->path()) {
        my $candidate = File::Spec->catfile($dir, 'starfix_cli');
        return abs_path($candidate) if -f $candidate && -x $candidate;
    }
    die "starfix_cli was not found; run `make starfix_cli`\n";
}

sub get_db_path {
    for my $candidate ('data/starfix_db.bin', '../data/starfix_db.bin') {
        if (-f $candidate) {
            return abs_path($candidate);
        }
    }
    die "data/starfix_db.bin was not found\n";
}

sub estimate_crop_fov {
    my ($exif_dict, $width, $height, $crop_size) = @_;
    $crop_size //= $CROP_SIZE;
    my $exif_ifd = ($exif_dict && $exif_dict->{'Exif'}) ? $exif_dict->{'Exif'} : {};
    my $focal_35 = $exif_ifd->{$TAG_FOCAL_LENGTH_35MM};
    if (!defined $focal_35 || !is_finite_number($focal_35) || $focal_35 <= 0.0) {
        $focal_35 = 24.0;
    } else {
        $focal_35 += 0.0;
    }

    # 35 mm equivalent focal lengths preserve diagonal FOV, including on 4:3 sensors
    my $diag_pixels = sqrt($width * $width + $height * $height);
    my $diag_35mm = sqrt(36.0 * 36.0 + 24.0 * 24.0);
    my $focal_pixels = $diag_pixels * $focal_35 / $diag_35mm;
    return 2.0 * (atan2($crop_size / 2.0, $focal_pixels) * 180.0 / $PI);
}

sub run_command {
    my ($command_ref, $work_dir) = @_;
    my $pid = open(my $fh, '-|');
    die "fork failed: $!\n" unless defined $pid;
    if ($pid == 0) {
        if (defined $work_dir) {
            chdir($work_dir) or die "chdir $work_dir: $!\n";
        }
        open(STDERR, '>&', STDOUT) or die "redirect stderr: $!\n";
        exec { $command_ref->[0] } @$command_ref;
        die "exec failed: $!\n";
    }
    local $/;
    my $stdout = <$fh>;
    close($fh);
    my $status = $?;
    my $returncode = ($status == -1) ? -1 : ($status >> 8);
    return {
        stdout     => defined $stdout ? $stdout : '',
        stderr     => '',
        returncode => $returncode,
    };
}

sub parse_solver_output {
    my ($result) = @_;
    my %parsed = (
        rms       => 9**9**9,
        rotation  => undef,
        inliers   => 0,
        matches   => 0,
        fov_deg   => undef,
        latitude  => undef,
        longitude => undef,
        stdout    => defined($result) ? ($result->{stdout} // '') : '',
    );
    return \%parsed if !defined $result || ($result->{returncode} // 1) != 0;

    for my $line (split(/\r?\n/, $result->{stdout} // '')) {
        if ($line =~ /\ASOLVED_FOV:\s+(\S+)/) {
            $parsed{fov_deg} = $1 + 0.0 if is_finite_number($1);
        } elsif ($line =~ /\ARMS_ERROR:\s+(\S+)/) {
            $parsed{rms} = $1 + 0.0 if is_finite_number($1);
        } elsif ($line =~ /\AINLIERS:\s+(\d+)\/(\d+)\z/) {
            $parsed{inliers} = int($1);
            $parsed{matches} = int($2);
        } elsif ($line =~ /\AROTATION_MATRIX:\s+(.+)\z/) {
            my @parts = split(/\s+/, $1);
            if (@parts == 9 && (grep { is_finite_number($_) } @parts) == 9) {
                my @v = map { $_ + 0.0 } @parts;
                $parsed{rotation} = [
                    [ $v[0], $v[1], $v[2] ],
                    [ $v[3], $v[4], $v[5] ],
                    [ $v[6], $v[7], $v[8] ],
                ];
            }
        } elsif ($line =~ /Solved Latitude\s*:\s*([+-]?\d+(?:\.\d+)?)/) {
            $parsed{latitude} = $1 + 0.0;
        } elsif ($line =~ /Solved Longitude\s*:\s*([+-]?\d+(?:\.\d+)?)/) {
            $parsed{longitude} = $1 + 0.0;
        }
    }
    return \%parsed;
}

sub solution_is_acceptable {
    my ($solution) = @_;
    return 0 unless $solution && ref($solution) eq 'HASH';
    return 0 unless defined $solution->{rotation};
    return 0 unless is_finite_number($solution->{rms}) && $solution->{rms} <= $MAX_RMS_PIXELS;
    return 0 unless ($solution->{inliers} // 0) >= $MIN_INLIERS;
    return 0 unless ($solution->{matches} // 0) >= $MIN_INLIERS;
    return 0 unless ($solution->{inliers} / $solution->{matches}) >= 0.6;
    return 1;
}

sub solve_plate {
    my ($image_path, $exif_dict, $cli_bin, $work_dir, $utc_dt, $zenith_camera) = @_;
    $zenith_camera ||= [ 0.0, 0.0, 1.0 ];
    my ($width, $height) = get_oriented_dimensions($image_path, $exif_dict);
    my $min_dim = $width < $height ? $width : $height;
    if ($min_dim < $CROP_SIZE) {
        die "image is smaller than the required ${CROP_SIZE}x${CROP_SIZE} crop\n";
    }
    my $year = decimal_year($utc_dt);
    my (undef, $gha_aries) = calculate_gha_aries($utc_dt);
    my $data_dir = File::Spec->catdir($work_dir, 'data');
    my $matches_path = File::Spec->catfile($data_dir, 'identified_stars.json');
    my $centroids_path = File::Spec->catfile($data_dir, 'detected_centroids.json');
    my $pgm_path = File::Spec->catfile($data_dir, 'mock_image.pgm');

    my %size_set;
    for my $factor (1, 1.5, 2, 2.5, 3, 4) {
        my $candidate = int($CROP_SIZE * $factor);
        $candidate = $min_dim if $candidate > $min_dim;
        $size_set{$candidate} = 1;
    }
    my @sizes = sort { $a <=> $b } keys %size_set;

    for my $crop_size (@sizes) {
        crop_center_to_pgm($image_path, $pgm_path, $crop_size);
        my $fov_estimate = estimate_crop_fov($exif_dict, $width, $height, $crop_size);

        for my $threshold (45, 35) {
            my $detected_res = run_command([ $cli_bin, '--centroid', "$threshold" ], $work_dir);
            next if $detected_res->{returncode} != 0;

            open(my $cfh, '<:encoding(UTF-8)', $centroids_path) or next;
            local $/;
            my $detected = decode_json(<$cfh>);
            close($cfh);

            # reject edge fragments and low-flux noise before selecting brightest quads
            my @centroids = grep {
                $_->{flux} >= 100
                    && $_->{u} >= 15 && $_->{u} <= $crop_size - 15
                    && $_->{v} >= 15 && $_->{v} <= $crop_size - 15
            } @$detected;
            printf("  Crop %dpx, threshold %d: %d usable stars\n",
                   $crop_size, $threshold, scalar(@centroids));
            next if @centroids < 5;

            open(my $ofh, '>:encoding(UTF-8)', $centroids_path)
                or die "cannot write $centroids_path: $!\n";
            print $ofh encode_json(\@centroids);
            close($ofh);

            my $sweep_res = run_command([
                $cli_bin, '--sweep-photo',
                sprintf('%.9f', $gha_aries),
                sprintf('%.12f', $zenith_camera->[0]),
                sprintf('%.12f', $zenith_camera->[1]),
                sprintf('%.12f', $zenith_camera->[2]),
                sprintf('%.9f', $fov_estimate),
                sprintf('%.9f', $year),
                "$crop_size",
            ], $work_dir);
            my $solution = parse_solver_output($sweep_res);
            if (solution_is_acceptable($solution) && defined $solution->{fov_deg}) {
                open(my $mfh, '<:encoding(UTF-8)', $matches_path)
                    or die "cannot read $matches_path: $!\n";
                local $/;
                my $identified = decode_json(<$mfh>);
                close($mfh);
                return {
                    crop_size  => $crop_size,
                    fov_deg    => $solution->{fov_deg},
                    solution   => $solution,
                    centroids  => \@centroids,
                    identified => $identified,
                    width      => $width,
                    height     => $height,
                };
            }
        }
    }
    die "no verified plate solution after searching central and wider crops\n";
}

1;
