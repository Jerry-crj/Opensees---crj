proc divideStepAnalysis {DT  divideNum  tol} {
	puts "function is called..."
	set  minDt  5.e-6
    set  dt  [ expr   $DT/$divideNum ]
	if {abs($dt) > abs($minDt)} {
		for {set i 0} {$i < $divideNum} {incr i} {
			set  newTol  [expr $tol]
			test NormDispIncr  $newTol  10  1
			integrator LoadControl $dt
			puts "dt = $dt, time = [getTime], divide step analysis"
			set ok [analyze 1]
			if {$ok != 0} {
				set ok [divideStepAnalysis $dt $divideNum  $newTol]
			}
		}
	} else {
		integrator LoadControl $DT
		puts "dt = $DT, time = [getTime], using other algorithm ... "
		puts "algorithm NewtonLineSearch"
		algorithm NewtonLineSearch 
		set ok [analyze 1] 
		if {$ok != 0} {				
			puts "algorithm inital Newton"
			test NormDispIncr $tol  20  5
			algorithm Newton -initial
			analyze 1
		}
	}
	return $ok
}