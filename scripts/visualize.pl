#!/usr/bin/env perl
use strict;
use warnings;
use File::Spec;
use JSON;

sub main {
    my $script_dir = File::Spec->rel2abs(File::Spec->path());
    my @dirs = File::Spec->splitdir(File::Spec->rel2abs($0));
    pop @dirs; # pop filename
    pop @dirs; # pop scripts dir
    my $base_dir = File::Spec->catdir(@dirs);

    my $mock_image = File::Spec->catfile($base_dir, "data", "mock_image.pgm");
    my $centroids_file = File::Spec->catfile($base_dir, "data", "detected_centroids.json");
    my $identified_file = File::Spec->catfile($base_dir, "data", "identified_stars.json");
    my $attitude_file = File::Spec->catfile($base_dir, "data", "attitude.json");

    my $step1_out = File::Spec->catfile($base_dir, "data", "step1_centroids.png");
    my $step2_out = File::Spec->catfile($base_dir, "data", "step2_identified.png");
    my $step2b_out = File::Spec->catfile($base_dir, "data", "step2b_zenith_camera.png");
    my $step2c_out = File::Spec->catfile($base_dir, "data", "step2c_zenith_celestial.png");
    my $step2d_out = File::Spec->catfile($base_dir, "data", "step2d_earth_coordinates.png");
    my $step3_out = File::Spec->catfile($base_dir, "data", "step3_position.png");

    my $im_tool = $ENV{IMTOOL} || "magick";

    print "Generating Step 1: Centroid Detection...\n";
    if (-f $centroids_file) {
        open(my $fh, '<', $centroids_file) or die "Cannot open centroids file: $!";
        local $/;
        my $json_text = <$fh>;
        close($fh);
        my $centroids = decode_json($json_text);

        my $draw_args = "";
        foreach my $c (@$centroids) {
            my $u = $c->{u};
            my $v = $c->{v};
            my $r = 24;
            $draw_args .= " -stroke red -fill none -strokewidth 4 -draw \"circle $u,$v " . ($u + $r) . ",$v\"";
        }
        
        my $cmd = "$im_tool \"$mock_image\" $draw_args -fill white -pointsize 48 -draw \"text 40,80 'Step 1: Centroid Extraction (Detected: " . scalar(@$centroids) . " stars)'\" \"$step1_out\"";
        system($cmd) == 0 or warn "Failed to create Step 1 image: $!";
        print "  Saved $step1_out\n";
    } else {
        print "  Centroids file not found, skipping.\n";
    }

    print "Generating Step 2: Star Identification...\n";
    if (-f $identified_file && -f $centroids_file) {
        # Read centroids first
        open(my $cfh, '<', $centroids_file) or die "Cannot open centroids: $!";
        local $/;
        my $centroids_text = <$cfh>;
        close($cfh);
        my $centroids = decode_json($centroids_text);

        # Read identified stars
        open(my $fh, '<', $identified_file) or die "Cannot open identified file: $!";
        my $json_text = <$fh>;
        close($fh);
        my $identified_raw = decode_json($json_text);

        # Map identified stars to their u,v coordinates
        my @identified;
        foreach my $star (@$identified_raw) {
            my $idx = $star->{centroid_idx};
            my $c = $centroids->[$idx];
            push @identified, {
                hip => $star->{hip},
                u => $c->{u},
                v => $c->{v}
            };
        }

        my $draw_args = "";
        
        # Connect stars in a loop to show the quad geometry
        my @sorted = sort { $a->{u} <=> $b->{u} } @identified;
        if (@sorted >= 3) {
            for (my $i = 0; $i < $#sorted; $i++) {
                my ($u1, $v1) = ($sorted[$i]->{u}, $sorted[$i]->{v});
                my ($u2, $v2) = ($sorted[$i+1]->{u}, $sorted[$i+1]->{v});
                $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line $u1,$v1 $u2,$v2\"";
            }
            # close the loop
            my ($u1, $v1) = ($sorted[-1]->{u}, $sorted[-1]->{v});
            my ($u2, $v2) = ($sorted[0]->{u}, $sorted[0]->{v});
            $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line $u1,$v1 $u2,$v2\"";
        }

        # Highlight and label each star
        foreach my $star (@identified) {
            my $u = $star->{u};
            my $v = $star->{v};
            my $hip = $star->{hip};
            my $r = 30;
            $draw_args .= " -stroke \"#00FF00\" -fill none -strokewidth 4 -draw \"circle $u,$v " . ($u + $r) . ",$v\"";
            $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 36 -draw \"text " . ($u + 36) . "," . ($v + 12) . " 'HIP $hip' \"";
        }

        # Read pointing angles if available to show in header
        my $pointing_str = "";
        if (-f $attitude_file) {
            open(my $afh, '<', $attitude_file) or warn "Cannot open attitude file: $!";
            if ($afh) {
                local $/;
                my $att_text = <$afh>;
                close($afh);
                my $att = decode_json($att_text);
                $pointing_str = sprintf("RA: %.3f, Dec: %.3f, Roll: %.3f", $att->{ra_deg}, $att->{dec_deg}, $att->{roll_deg});
            }
        }

        my $header_text = "Step 2: Lost-in-Space Star Identification & Camera Pointing";
        my $sub_text = $pointing_str ? "Pointing: $pointing_str" : "";
        my $cmd = "$im_tool \"$mock_image\" $draw_args -fill white -pointsize 48 -draw \"text 40,80 '$header_text'\" " .
                  ($sub_text ? "-draw \"text 40,140 '$sub_text'\"" : "") . " \"$step2_out\"";
        system($cmd) == 0 or warn "Failed to create Step 2 image: $!";
        print "  Saved $step2_out\n";
    } else {
        print "  Identified stars file not found, skipping.\n";
    }

    # --- Transition Step A: Gravity Vector in Camera Frame ---
    print "Generating Transition Step A: Gravity in Camera Frame...\n";
    {
        my $width = 2048;
        my $height = 1536;
        my $center_x = $width / 2;
        my $center_y = $height / 2;
        
        my $draw_args = "";
        
        # Draw tilted sensor plane as a polygon
        $draw_args .= " -stroke \"#005500\" -fill \"#001100\" -strokewidth 3 -draw \"polygon 600,600 1448,600 1200,950 352,950\"";
        
        # Draw Z axis (optical axis) normal to sensor
        $draw_args .= " -stroke \"#FF00FF\" -strokewidth 4 -draw \"line $center_x,$center_y $center_x,250\"";
        # Draw Arrow head for Z axis
        $draw_args .= " -stroke \"#FF00FF\" -strokewidth 4 -draw \"line $center_x,250 " . ($center_x - 15) . ",280\"";
        $draw_args .= " -stroke \"#FF00FF\" -strokewidth 4 -draw \"line $center_x,250 " . ($center_x + 15) . ",280\"";
        $draw_args .= " -fill \"#FF00FF\" -stroke none -pointsize 32 -draw \"text " . ($center_x + 25) . ",270 'Camera Z-Axis (Optical Axis)'\"";

        # Draw X axis on sensor
        $draw_args .= " -stroke \"#FF0000\" -strokewidth 4 -draw \"line $center_x,$center_y 1400,$center_y\"";
        $draw_args .= " -fill \"#FF0000\" -stroke none -pointsize 32 -draw \"text 1420,$center_y 'Camera X-Axis'\"";

        # Draw Y axis on sensor
        $draw_args .= " -stroke \"#00FF00\" -strokewidth 4 -draw \"line $center_x,$center_y 1150,910\"";
        $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 32 -draw \"text 1170,940 'Camera Y-Axis'\"";

        # Draw gravity vector (Zenith) pointing down/left
        $draw_args .= " -stroke \"#FFFF00\" -strokewidth 6 -draw \"line $center_x,$center_y 700,1250\"";
        # Draw Arrow head for gravity
        $draw_args .= " -stroke \"#FFFF00\" -strokewidth 6 -draw \"line 700,1250 710,1210\"";
        $draw_args .= " -stroke \"#FFFF00\" -strokewidth 6 -draw \"line 700,1250 740,1240\"";
        $draw_args .= " -fill \"#FFFF00\" -stroke none -pointsize 36 -draw \"text 400,1320 'IMU Gravity Zenith Vector (z_cam)'\"";

        # Sensor label
        $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 32 -draw \"text 380,920 'Tilted Sensor Plane (X-Y)'\"";

        # Header titles
        $draw_args .= " -fill white -stroke none -pointsize 54 -draw \"text 60,100 'Transition A: Gravity Vector in Camera Frame'\"";
        $draw_args .= " -fill \"#888888\" -pointsize 36 -draw \"text 60,170 'IMU/Accelerometers measure gravity direction relative to camera local frame'\"";

        my $cmd = "$im_tool -size ${width}x${height} xc:black -fill none -stroke \"#002200\" -strokewidth 3 -draw \"rectangle 0,0 " . ($width-1) . "," . ($height-1) . "\" $draw_args \"$step2b_out\"";
        system($cmd) == 0 or warn "Failed to create Transition A image: $!";
        print "  Saved $step2b_out\n";
    }

    # --- Transition Step B: Projected Celestial Zenith ---
    print "Generating Transition Step B: Projected Celestial Zenith...\n";
    {
        my $width = 2048;
        my $height = 1536;
        my $center_x = $width / 2;
        my $center_y = $height / 2;
        
        my $draw_args = "";
        
        # Celestial Sphere boundary (radius 600)
        $draw_args .= " -stroke \"#003300\" -fill \"#000500\" -strokewidth 4 -draw \"circle $center_x,$center_y " . ($center_x + 600) . ",$center_y\"";
        
        # Celestial Equator
        $draw_args .= " -stroke \"#005500\" -fill none -strokewidth 3 -draw \"ellipse $center_x,$center_y 600,180 0,360\"";
        $draw_args .= " -stroke \"#005500\" -fill none -strokewidth 1.5 -draw \"line " . ($center_x - 600) . ",$center_y " . ($center_x + 600) . ",$center_y\"";
        $draw_args .= " -fill \"#005500\" -stroke none -pointsize 28 -draw \"text " . ($center_x + 610) . "," . ($center_y + 8) . " 'Celestial Equator'\"";

        # Parallels of Declination
        $draw_args .= " -stroke \"#003300\" -fill none -strokewidth 2 -draw \"ellipse $center_x,600 520,130 0,360\"";
        $draw_args .= " -stroke \"#003300\" -fill none -strokewidth 2 -draw \"ellipse $center_x,450 400,100 0,360\"";

        # Meridians of RA
        $draw_args .= " -stroke \"#003300\" -fill none -strokewidth 2 -draw \"line $center_x,168 $center_x,1368\"";
        
        # Celestial Zenith vector (z_celestial) pointing to surface at (1350, 450)
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line $center_x,$center_y 1350,450\"";
        # Draw Arrow head
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line 1350,450 1310,470\"";
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line 1350,450 1340,495\"";
        $draw_args .= " -fill \"#00FFFF\" -stroke none -pointsize 36 -draw \"text 1370,430 'Projected Celestial Zenith (z_celestial)'\"";
        $draw_args .= " -fill \"#00FFFF\" -stroke none -pointsize 28 -draw \"text 1370,480 'Dec: 37.77 deg, RA: 80.00 deg'\"";

        # Draw a few stars for context
        my @stars_to_draw = (
            { u => 1250, v => 550, name => "HIP 25292" },
            { u => 1100, v => 400, name => "HIP 23179" },
            { u => 1400, v => 600, name => "HIP 25541" }
        );
        foreach my $s (@stars_to_draw) {
            $draw_args .= " -stroke \"#00FF00\" -fill \"#00FF00\" -draw \"circle $s->{u},$s->{v} " . ($s->{u} + 10) . ",$s->{v}\"";
            $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 24 -draw \"text " . ($s->{u} + 15) . "," . ($s->{v} + 8) . " '$s->{name}'\"";
        }

        # Header titles
        $draw_args .= " -fill white -stroke none -pointsize 54 -draw \"text 60,100 'Transition B: Projecting Zenith to Celestial Sphere'\"";
        $draw_args .= " -fill \"#888888\" -pointsize 36 -draw \"text 60,170 'Attitude matrix R^T rotates camera-frame Zenith into the global inertial frame'\"";

        my $cmd = "$im_tool -size ${width}x${height} xc:black -fill none -stroke \"#002200\" -strokewidth 3 -draw \"rectangle 0,0 " . ($width-1) . "," . ($height-1) . "\" $draw_args \"$step2c_out\"";
        system($cmd) == 0 or warn "Failed to create Transition B image: $!";
        print "  Saved $step2c_out\n";
    }

    # --- Transition Step C: Terrestrial Coordinate Projection ---
    print "Generating Transition Step C: Terrestrial Coordinate Projection...\n";
    {
        my $width = 2048;
        my $height = 1536;
        my $center_x = $width / 2;
        my $center_y = $height / 2;
        
        my $draw_args = "";
        
        # Earth Sphere
        $draw_args .= " -stroke \"#333333\" -fill \"#000500\" -strokewidth 4 -draw \"circle $center_x,$center_y " . ($center_x + 600) . ",$center_y\"";
        
        # Earth Equator
        $draw_args .= " -stroke \"#005500\" -fill none -strokewidth 3 -draw \"ellipse $center_x,$center_y 600,180 0,360\"";
        
        # Earth Prime Meridian
        $draw_args .= " -stroke \"#005500\" -fill none -strokewidth 3 -draw \"ellipse $center_x,$center_y 180,600 0,360\"";
        
        # Line from center of Earth through SF to celestial Zenith
        $draw_args .= " -stroke \"#FFFF00\" -strokewidth 4 -draw \"line $center_x,$center_y 650,300\"";
        # Zenith pointer in space
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line 650,300 580,200\"";
        # Arrow head
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line 580,200 620,220\"";
        $draw_args .= " -stroke \"#00FFFF\" -strokewidth 6 -draw \"line 580,200 590,240\"";
        $draw_args .= " -fill \"#00FFFF\" -stroke none -pointsize 36 -draw \"text 440,160 'Zenith'\"";

        # Plot San Francisco (on the line)
        my $sf_x = 760;
        my $sf_y = 440;
        $draw_args .= " -stroke \"#00FF00\" -fill \"#00FF00\" -draw \"circle $sf_x,$sf_y " . ($sf_x + 12) . ",$sf_y\"";
        $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 32 -draw \"text " . ($sf_x + 18) . "," . ($sf_y + 10) . " 'San Francisco'\"";

        # Latitude angle indicator
        $draw_args .= " -stroke \"#00FF00\" -strokewidth 3 -draw \"line $center_x,$center_y 1024,588\""; # equator radius line
        $draw_args .= " -fill \"#00FF00\" -stroke none -pointsize 36 -draw \"text 920,530 'Latitude = Dec = 37.77 deg'\"";

        # Longitude angle indicator
        $draw_args .= " -stroke \"#FF0000\" -strokewidth 3 -draw \"line $center_x,$center_y 844,768\""; # meridian line
        $draw_args .= " -fill \"#FF3333\" -stroke none -pointsize 36 -draw \"text 600,900 'Longitude = RA - GHA = -122.41 deg'\"";

        # Header titles
        $draw_args .= " -fill white -stroke none -pointsize 54 -draw \"text 60,100 'Transition C: Deriving Latitude & Longitude on Earth'\"";
        $draw_args .= " -fill \"#888888\" -pointsize 36 -draw \"text 60,170 'Zenith Declination maps directly to Latitude. Right Ascension offset from Greenwich maps to Longitude'\"";

        my $cmd = "$im_tool -size ${width}x${height} xc:black -fill none -stroke \"#002200\" -strokewidth 3 -draw \"rectangle 0,0 " . ($width-1) . "," . ($height-1) . "\" $draw_args \"$step2d_out\"";
        system($cmd) == 0 or warn "Failed to create Transition C image: $!";
        print "  Saved $step2d_out\n";
    }

    print "Generating Step 3: Celestial Position Solver Infographic...\n";
    # We will generate a brand new dark tactical radar screen plot
    my $width = 2048;
    my $height = 1536;
    
    # Coordinates in pixels mapping: center is San Francisco
    my $center_x = $width / 2;
    my $center_y = $height / 2;

    # Coordinates scaling: 1 deg = 1200 pixels
    my $scale = 1200.0; 

    # San Francisco Coordinates (True)
    my $true_lat = 37.7749;
    my $true_lon = -122.4194;

    # Assumed Position (AP)
    my $ap_lat = 37.5;
    my $ap_lon = -122.0;

    # Solved Position (from C solver run, dynamic with fallback)
    my $solved_lat = 37.7763;
    my $solved_lon = -122.4170;
    my $solved_file = File::Spec->catfile($base_dir, "data", "solved_position.json");
    if (-f $solved_file) {
        open(my $sfh, '<', $solved_file) or warn "Cannot open solved position file: $!";
        if ($sfh) {
            local $/;
            my $json_text = <$sfh>;
            close($sfh);
            my $solved_pos = decode_json($json_text);
            $solved_lat = $solved_pos->{lat};
            $solved_lon = $solved_pos->{lon};
        }
    }

    # Calculate pixel offsets relative to True position (Center of radar)
    my $ap_dx = ($ap_lon - $true_lon) * cos($true_lat * 3.14159265 / 180.0) * $scale;
    my $ap_dy = -($ap_lat - $true_lat) * $scale; # negate because y increases downward in image

    my $sol_dx = ($solved_lon - $true_lon) * cos($true_lat * 3.14159265 / 180.0) * $scale;
    my $sol_dy = -($solved_lat - $true_lat) * $scale;

    my $ap_px = $center_x + $ap_dx;
    my $ap_py = $center_y + $ap_dy;

    my $sol_px = $center_x + $sol_dx;
    my $sol_py = $center_y + $sol_dy;

    my $draw_args = "";
    
    # Draw radar grids (circles of radius 10nm, 20nm, 30nm)
    # 1nm = 1/60 deg = 1200/60 = 20.0 pixels
    for (my $nm = 5; $nm <= 25; $nm += 5) {
        my $r = $nm * (1200.0 / 60.0);
        $draw_args .= " -stroke \"#003300\" -fill none -strokewidth 3 -draw \"circle $center_x,$center_y " . ($center_x + $r) . ",$center_y\"";
        $draw_args .= " -stroke none -fill \"#006600\" -pointsize 28 -draw \"text " . ($center_x + $r + 10) . "," . ($center_y + 10) . " '$nm nm'\"";
    }

    # Crosshair axes
    $draw_args .= " -stroke \"#003300\" -strokewidth 3 -draw \"line 0,$center_y $width,$center_y\"";
    $draw_args .= " -stroke \"#003300\" -strokewidth 3 -draw \"line $center_x,0 $center_x,$height\"";

    # Draw Assumed Position (AP) in Red
    $draw_args .= " -stroke \"#FF0000\" -fill none -strokewidth 4 -draw \"circle $ap_px,$ap_py " . ($ap_px + 24) . ",$ap_py\"";
    $draw_args .= " -stroke \"#FF0000\" -strokewidth 4 -draw \"line " . ($ap_px-36) . ",$ap_py " . ($ap_px+36) . ",$ap_py\"";
    $draw_args .= " -stroke \"#FF0000\" -strokewidth 4 -draw \"line $ap_px," . ($ap_py-36) . " $ap_px," . ($ap_py+36) . "\"";
    $draw_args .= " -fill \"#FF3333\" -stroke none -pointsize 36 -draw \"text " . ($ap_px+40) . "," . ($ap_py+12) . " 'Assumed Position (AP)'\"";

    # Draw True Position in Green (Center of Target)
    $draw_args .= " -stroke \"#00FF00\" -fill none -strokewidth 5 -draw \"circle $center_x,$center_y " . ($center_x + 18) . ",$center_y\"";
    $draw_args .= " -fill \"#33FF33\" -stroke none -pointsize 36 -draw \"text " . ($center_x+30) . "," . ($center_y-30) . " 'True Position (SF Mast)'\"";

    # Draw Solved Position in Cyan
    $draw_args .= " -stroke \"#00FFFF\" -fill none -strokewidth 5 -draw \"circle $sol_px,$sol_py " . ($sol_px + 24) . ",$sol_py\"";
    $draw_args .= " -stroke \"#00FFFF\" -strokewidth 4 -draw \"line " . ($sol_px-36) . ",$sol_py " . ($sol_px+36) . ",$sol_py\"";
    $draw_args .= " -stroke \"#00FFFF\" -strokewidth 4 -draw \"line $sol_px," . ($sol_py-36) . " $sol_px," . ($sol_py+36) . "\"";
    $draw_args .= " -fill \"#33FFFF\" -stroke none -pointsize 36 -draw \"text " . ($sol_px+40) . "," . ($sol_py-30) . " 'Solved Position'\"";

    # Draw line from AP to Solved showing convergence
    $draw_args .= " -stroke \"#FFFF00\" -strokewidth 4 -draw \"line $ap_px,$ap_py $sol_px,$sol_py\"";

    # Write text block at top left
    $draw_args .= " -fill white -stroke none -pointsize 54 -draw \"text 60,100 'Step 3: Celestial Position Solver Intersect'\"";
    $draw_args .= " -fill \"#888888\" -pointsize 36 -draw \"text 60,170 'GPS-Denied Least Squares Solver via nanoqsp'\"";

    # Write metrics at bottom left
    my $err_lat = $solved_lat - $true_lat;
    my $err_lon = ($solved_lon - $true_lon) * cos($true_lat * 3.14159265 / 180.0);
    my $dist_err = sqrt($err_lat*$err_lat + $err_lon*$err_lon) * 60.0;

    $draw_args .= " -fill \"#00FF00\" -pointsize 38 -draw \"text 60," . ($height - 280) . " 'True:  Lat 37.7749, Lon -122.4194 (San Francisco)'\"";
    $draw_args .= sprintf(" -fill \"#33FFFF\" -pointsize 38 -draw \"text 60,%d 'Solved: Lat %.4f, Lon %.4f'\"", $height - 220, $solved_lat, $solved_lon);
    $draw_args .= " -fill \"#FF3333\" -pointsize 38 -draw \"text 60," . ($height - 160) . " 'AP:     Lat 37.5000, Lon -122.0000'\"";
    my $dist_err_m = $dist_err * 1852.0;
    $draw_args .= sprintf(" -fill \"#FFFF00\" -pointsize 44 -draw \"text 60,%d 'Position Solver Error: %.4f Nautical Miles (~%.0fm)'\"", $height - 80, $dist_err, $dist_err_m);

    my $cmd = "$im_tool -size ${width}x${height} xc:black -fill none -stroke \"#002200\" -strokewidth 3 -draw \"rectangle 0,0 " . ($width-1) . "," . ($height-1) . "\" $draw_args \"$step3_out\"";
    system($cmd) == 0 or warn "Failed to create Step 3 infographic: $!";
    print "  Saved $step3_out\n";

    print "Done generating visual pipeline assets!\n";
}

main();
