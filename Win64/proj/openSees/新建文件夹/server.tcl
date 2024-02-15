 
# proc analyzeModel {data} {
	# # puts "[llength $data]"
	# # puts "$data"
	
	# # for {set i 0} {$i < [llength $data]} {incr i} {
		# # puts "$data"
		
		# source load.tcl

		# set ok [analyze 1]
		# if {$ok != 0} {
			# divideStepAnalysis   $dt    10    $tol
		# }
		# global times
		# global Fx Fy Fz
		# incr times
		# source RnodeTag.tcl
		# source LnodeTag.tcl

		# global F
		
		# for {set i 0} {$i < [llength $RnodeTag]} {incr i} {
			# set tagR [lindex $RnodeTag $i]
			# set tagL [lindex $LnodeTag $i]
			# set Fx [expr ([nodeDisp $tagL 1]*1e10+[nodeReaction $tagR 1])]
			# set Fy [expr ([nodeDisp $tagL 2]*1e10+[nodeReaction $tagR 2])]
			# set Fz [expr ([nodeDisp $tagL 3]*1e10+[nodeReaction $tagR 3])]
			# lappend  F $Fx $Fy $Fz
		# }
		
		# puts "$F"
		
		# puts "!!!!!!!!!!!!!!!!!!!!!"
		# # puts "Fx = $Fx; Fy = $Fy; Fz = $Fz"
	
		# # return
		# # source recordEachFrameNodeDispl0.tcl
	# # }
# }
 
proc accept {sock ip port} {
	fconfigure $sock -blocking 1 -buffering none ;#line
	fileevent $sock readable [list respond $sock]
}
 
proc respond {sock} {
	if {[eof $sock] || [catch {gets $sock data}]} {
		# end of file or abnormal connection drop
		close $sock
		puts "Close $sock" ;# $echo(addr,$sock)"
	} else {
		# safe calculator style evaluator
		set F {}
		source load.tcl

		set ok [analyze 1]
		if {$ok != 0} {
			divideStepAnalysis   $dt    10    $tol
		}
		global times
		# global Fx Fy Fz
		incr times
		source RnodeTag.tcl
		source LnodeTag.tcl

		# global F
		
		# for {set i 0} {$i < [llength $RnodeTag]} {incr i} {
			# set tagR [lindex $RnodeTag $i]
			# set tagL [lindex $LnodeTag $i]
			# set Fx [expr ([nodeDisp $tagL 1]*1e10+[nodeReaction $tagR 1])]
			# set Fy [expr ([nodeDisp $tagL 2]*1e10+[nodeReaction $tagR 2])]
			# set Fz [expr ([nodeDisp $tagL 3]*1e10+[nodeReaction $tagR 3])]
			# lappend  F $Fx $Fy $Fz
		# }

		for {set i 0} {$i < [llength $RnodeTag]} {incr i} {
			set tagR [lindex $RnodeTag $i]
			set tagL [lindex $LnodeTag $i]
			set Fx [nodeReaction $tagR 1]
			set Fy [nodeReaction $tagR 2]
			set Fz [nodeReaction $tagR 3]
			puts  "Reaction  $tagR  $Fx  $Fy  $Fz"
			lappend  F $Fx $Fy $Fz
		}

		
		# puts  "$F"
		puts $sock  "$F"
		flush $sock
		return
	}
}

# ---- model ----------

source mesosub.tcl
socket -server accept 7271
vwait forever