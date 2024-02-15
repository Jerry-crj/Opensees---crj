
source LnodeTag.tcl

remove loadPattern 1

# set u {}
pattern  Plain  1  1 {
	for {set i 0} {$i < [llength $LnodeTag]} {incr i} {
		set  tag   [lindex $LnodeTag $i]
		# set   ux   [lindex $data [expr 3*$i+0]]
		# set   uy   [lindex $data [expr 3*$i+1]]
		# set   uz   [lindex $data [expr 3*$i+2]]
		puts  "load  $tag    $ux    $uy    $uz"
		load  $tag    $ux    $uy    $uz
		# lappend    u    $ux    $uy    $uz
	}
}
# puts $u
