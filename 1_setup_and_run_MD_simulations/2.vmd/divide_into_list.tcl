lappend auto_path Orient
lappend auto_path la1.0

package require Orient
namespace import Orient::orient

set script  [open script w]

exec rm serials serials_tmp serials_tmp1 serials_ok
puts $script "sed -e 's/{//g' -e 's/}//g' -e 's/ /,/g' serials >> serials_ok"
#puts $script "sed 's/}//g' serials_tmp >> serials_tmp1"
#puts $script "sed 's/ /,/g' serials_tmp1 >> serials_ok"

close $script

exec chmod +x script
############################################  CHANGE HERE THE LIST OF RESIDUES LINING THE BS OF INTEREST. THE FOLLOWING ARE FOR ADK ########################
set sel [atomselect top "noh and resid 8 to 13 15 31 32 35 36 57 to 59 64 84 to 86 88 92 119 120 123 134 137 156 158 167 200 201 202 205"]
set selrog [atomselect top "[$sel text] and backbone"]
set L [lsort -integer -uniq [$sel get resid]]
set Rog_serials [lsort -integer -uniq [$selrog get serial]]
#print RoG serial numbers
set serials [open serials w]
puts $serials "Rog serials:"
puts $serials $Rog_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
mv serials_ok serials_RoG

############################################################################################################################################################
set cm_tmp [measure center $sel]
set cm_bs [vecscale -1 $cm_tmp]
set all [atomselect top all]
$all moveby $cm_bs
##ora allinea assi
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 2] {0 0 1}]
$all move $A
set I [draw principalaxes $sel]
set A [orient $sel [lindex $I 1] {0 1 0}]
$all move $A
set I [draw principalaxes $sel]

#plane xy
set Z_1  {}
set Z_2  {}
#plane xz
set Y_1  {}
set Y_2  {}
#plane yz
set X_1  {}
set X_2  {}

set Z_1_serials  {}
set Z_2_serials {}

set Y_1_serials  {}
set Y_2_serials {}

set X_1_serials  {}
set X_2_serials {}

foreach i $L {
    set resid [atomselect top "resid $i and backbone"]
    set cm [measure center $resid]

    if {[lindex $cm 0] >= 0} {
	lappend X_1 $i
    } else {lappend X_2 $i}
    if {[lindex $cm 1] >= 0} {
	lappend Y_1 $i
    } else {lappend Y_2 $i}
    if {[lindex $cm 2] >= 0} {
	lappend Z_1 $i
    } else {lappend Z_2 $i}
}

puts "plane xy"
puts "Z_1 resids: $Z_1"
puts "Z_2 resids: $Z_2"
puts ""
foreach i $Z_1 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend Z_1_serials [$resid_for_serial get serial]}
foreach i $Z_2 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend Z_2_serials [$resid_for_serial get serial]}

set serials [open serials w]
puts $serials "Z_1 serials:"
puts $serials $Z_1_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials

set serials [open serials w]
puts $serials "Z_2 serials:"
puts $serials $Z_2_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
exec mv serials_ok serials_xy

puts ""
puts "plane xz"
puts "Y_1 resids: $Y_1"
puts "Y_2 resids: $Y_2"

foreach i $Y_1 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend Y_1_serials [$resid_for_serial get serial]}
foreach i $Y_2 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend Y_2_serials [$resid_for_serial get serial]}
puts ""
set serials [open serials w]
puts $serials "Y_1 serials:"
puts $serials $Y_1_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
set serials [open serials w]
puts $serials "Y_2 serials:"
puts $serials $Y_2_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
exec mv serials_ok serials_xz

puts ""
puts "plane yz"
puts "X_1 resids: $X_1"
puts "X_2 resids: $X_2"

foreach i $X_1 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend X_1_serials [$resid_for_serial get serial]}
foreach i $X_2 {
    set resid_for_serial [atomselect top "resid $i and noh"]
    lappend X_2_serials [$resid_for_serial get serial]}
puts ""

set serials [open serials w]
puts $serials "X_1 serials:"
puts $serials $X_1_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
set serials [open serials w]
puts $serials "X_2 serials:"
puts $serials $X_2_serials
close $serials
exec ./script
puts [read [open serials_ok r]]
exec rm serials
exec mv serials_ok serials_yz

exec rm script

quit
