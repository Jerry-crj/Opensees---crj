
wipe

if { [file exists TimoFlexural] == 0 } {
	file mkdir TimoFlexural
}

logFile outputScreen.txt
model basic -ndm 2 -ndf 3

nDMaterial        ASIConcrete   1    -34     -1.7e-3     -5    -1.e-2    0.3    2.0   10.    2.0    0.006   1.0e-6
# uniaxialMaterial  Concrete02   103    -34     -1.7e-3     -5    -4.e-2    0.3    3.5   50.   
# uniaxialMaterial  Elastic      104    120
# nDMaterial        CombinedUniaxialTimoMat   1   103   104

uniaxialMaterial  Steel02      101    455      2.e5      0.001    14.0     0.925    0.15 
uniaxialMaterial  Hysteretic    102    28   3e-3   30  20e-3    30.01   22e-3   -28   -3e-3   -30  -20e-3    -30.01    -22e-3     0.0   0.0  0.0  0.0   0.4
   
nDMaterial        CombinedUniaxialTimoMat   2   101   102

node  1    0.0     0
node  2    0.0     137.083333333333
node  3    0.0     274.166666666667
node  4    0.0     411.25
node  5    0.0     548.333333333333
node  6    0.0     685.416666666667
node  7    0.0     822.5
node  8    0.0     959.583333333333
node  9    0.0     1096.66666666667
node  10   0.0     1233.75
node  11   0.0     1370.83333333333
node  12   0.0     1507.91666666667
node  13   0.0     1645

fix 1 1 1 1

set sectionNums 12
section FiberSectionTimoConstShear 1 {
	patch rect        1    $sectionNums      1      -175.    -175.     175.     175.
	
	layer straight    2    2     298.6477       -160.     -160.    -160.     160.    
	layer straight    2    2     298.6477        160.     -160.     160.     160.    
	layer straight    2    2     298.6477          0.     -160.       0.     160.
	layer straight    2    2     298.6477       -160.        0.     160.       0. 
}

element TimoFourNodePlaneSec2D  1  1  2  3  4  5  1 
element TimoFourNodePlaneSec2D  2  4  5  6  7  5  1 
element TimoFourNodePlaneSec2D  3  7  8  9  10  5  1 
element TimoFourNodePlaneSec2D  4  10 11 12 13  5  1 

uniaxialMaterial  Elastic  3   1.e14 

node  101    -1.0   1645.
fix 101 1 1 1

element truss  101  101  13  1.0  3

pattern Plain 1 Linear { 
    load  13  0.0  -1.0  0.0
}

system BandGeneral
test NormDispIncr 1.e-8 6
constraints Transformation
integrator LoadControl 83.1e3
algorithm Newton
numberer RCM
analysis Static
analyze 10

puts "gravity is finished..." 

# ------------------------------------------------------

loadConst -time 0 

recorder Node -file TimoFlexural/react.out    -precision 12 -time -node  1  -dof 1 2 3 reaction;
recorder Node -file TimoFlexural/displ.out    -precision 12 -time -node  13  -dof 1 2 3 disp;
# recorder Element -file TimoFlexural/fiber.out    -precision 12 -time -ele  1 2 3 4  section paperInfo;
recorder Element -file TimoFlexural/damage.out   -precision 12 -time -ele  1 2 3 4  section dmg

pattern Plain 2 Linear { 
    load  13  1.0e14   0.0   0.0
}

source GeneratePeaks.tcl
source divideStepAnalysis.tcl

set divideNum 10

# set iDmax "3"


set iDmax "15.39792388  -16.43598616	32.69896194	-33.39100346  48.61591696	-51.03806228	65.57093426	-67.64705882"
# HalfCycle

set tol 5.e-4

global iter;
set iter  1
foreach Dmax $iDmax {
	set DispCurve [GeneratePeaks  $Dmax  0.01  "HalfCycle"  1.0]
    set Displast 0.0
    foreach DispCurr $DispCurve {
	    set DispIncr [expr $DispCurr-$Displast] 
        system BandGeneral
        test NormDispIncr $tol 10
        constraints Transformation
        integrator LoadControl $DispIncr
        algorithm Newton 
        numberer RCM
        analysis Static
		set ok [analyze 1]
		if {$ok != 0} {
			divideStepAnalysis   $DispIncr  $divideNum  $tol
		}
		set Displast $DispCurr
	}
}

wipe



