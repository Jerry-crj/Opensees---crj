
wipe ;

set ux 0.0;
set uy 0.0;
set uz 1.0e2;

logFile screen.out
model basic -ndm 3 -ndf 3
set A 1.0
set times 0
set dt 1.0
set tol 3.0e-3

source material.tcl
source node.tcl
source fix.tcl
source element.tcl

timeSeries Constant 1 -factor 1.0e10
pattern  Plain  1  1 {
	load  3016  1.0   1.0   1.0e2
	# source load.tcl
}

# source reactRecorder0.tcl
# source displRecorder0.tcl
# source recordInitNodeCoord0.tcl
# source divideStepAnalysis.tcl

# constraints Transformation;     				
# numberer ParallelPlain
# system Mumps		
puts "1111"
constraints Transformation;     				
numberer Plain
system BandGeneral		
test NormDispIncr  $tol  10	 1	
algorithm Newton;					
integrator LoadControl  $dt
analysis Static	

puts "1111"

source RnodeTag.tcl
puts "1111"

for {set i 0} {$i < 5} {incr i} {
analyze 1
for {set i 0} {$i < [llength $RnodeTag]} {incr i} {
	set tagR [lindex $RnodeTag $i]
	# set tagL [lindex $LnodeTag $i]
	set Fx [nodeReaction $tagR 1]
	set Fy [nodeReaction $tagR 2]
	set Fz [nodeReaction $tagR 3]
	puts  "Reaction  $tagR  $Fx  $Fy  $Fz"
	lappend  F $Fx $Fy $Fz
}
}

# wipe