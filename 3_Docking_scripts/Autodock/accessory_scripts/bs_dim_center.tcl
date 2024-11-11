set out [open ./bs_dimension.dat w]
#puts $out "# frame        X              Y               Z                    center"
set gridspacing GGG
set num_steps [molinfo top get numframes]
set bs [atomselect top "XXX"]

for {set frame 0} {$frame < $num_steps} {incr frame} {
    $bs frame $frame
    $bs update
    set x_min [lindex [lindex [measure minmax $bs] 0] 0]
    set y_min [lindex [lindex [measure minmax $bs] 0] 1]
    set z_min [lindex [lindex [measure minmax $bs] 0] 2]
    set x_max [lindex [lindex [measure minmax $bs] 1] 0]
    set y_max [lindex [lindex [measure minmax $bs] 1] 1]
    set z_max [lindex [lindex [measure minmax $bs] 1] 2]
    
    set x_box [expr $x_max - $x_min]
    set y_box [expr $y_max - $y_min]
    set z_box [expr $z_max - $z_min]

    set nx [expr entier($x_box / $gridspacing)]
    set ny [expr entier($y_box / $gridspacing)]
    set nz [expr entier($z_box / $gridspacing)]
    
    set center [measure center $bs]
    set x_c [lindex $center 0]
    set y_c [lindex $center 1]
    set z_c [lindex $center 2]
    #puts $out "BS"

    #puts $out [format "%.3f %.3f %.3f %.3f %.3f %.3f" $x_box $y_box $z_box $x_c $y_c $z_c]
    puts $out [format "%d,%d,%d" $nx $ny $nz]
    puts $out [format "%.3f,%.3f,%.3f" $x_c $y_c $z_c]
}
close $out
quit
