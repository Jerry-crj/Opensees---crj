proc divideStepAnalysis {stepSize  divideNum  tol } {
	set  minStepSize  5.e-4
    set  newStepSize  [ expr   $stepSize/$divideNum ]
	if {abs($newStepSize) > abs($minStepSize)} {
		set  newTol   [expr $tol*$divideNum]
		test NormDispIncr $tol 10   1
		integrator LoadControl $newStepSize
		for {set i 0} {$i < $divideNum} {incr i} {
			set ok [analyze 1]
			if {$ok != 0} {
				set ok [divideStepAnalysis $newStepSize $divideNum  $tol]
			}
		}
	} else {
		test NormDispIncr $tol 1  5  
		analyze 1
		set ok 1
	}

	test NormDispIncr $tol 10  1
	integrator LoadControl $stepSize  $tol

	return $ok
}