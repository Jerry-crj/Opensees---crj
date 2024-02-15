
wipe

if { [file exists TimoShear] == 0 } {
	file mkdir TimoShear
}

logFile outputScreen.txt
model basic -ndm 2 -ndf 3

# uniaxialMaterial  Concrete02   103   -21.9    -3.0e-3     -2    -2.e-2    0.01    1.0   10. 
# uniaxialMaterial  Elastic      104    180
# nDMaterial        CombinedUniaxialTimoMat   1   103   104

nDMaterial        ASIConcrete   1    -21.9    -3.0e-3     -2    -2.e-2    0.01    0.5   1.     0.5    0.025   0.025 

uniaxialMaterial  Steel02      101     414      2.4e5  	    0.02    14.0     0.925    0.15 
uniaxialMaterial  Hysteretic    102    20   4e-3   25  20e-3    25.01   22e-3   -20   -4e-3   -25  -20e-3    -25.01    -22e-3     0.8  0.2  0.0  0.0   0.0
nDMaterial        CombinedUniaxialTimoMat   2   101   102

node   1     0.0    0
node   2     0.0    67.7333333333333
node   3     0.0    135.466666666667
node   4     0.0    203.2
node   5     0.0    270.933333333333
node   6     0.0    338.666666666667
node   7     0.0    406.4
node   8     0.0    474.133333333333
node   9     0.0    541.866666666667
node   10    0.0    609.6
node   11    0.0    677.333333333333
node   12    0.0    745.066666666667
node   13    0.0    812.8
node   14    0.0    880.533333333333
node   15    0.0    948.266666666667
node   16    0.0    1016
node   17    0.0    1083.73333333333
node   18    0.0    1151.46666666667
node   19    0.0    1219.2

fix 1 1 1 1

set sectionNums 12
section FiberSectionTimoConstShear  1   {
	patch rect       1    $sectionNums      1       -228.6   -457.2     228.6    457.2 
	
	layer straight   2    6    490.8739   -190.6     -419.2    -190.6     419.2
	layer straight   2    6    490.8739    190.6     -419.2     190.6     419.2
	layer straight   2    2    490.8739    -63.53    -419.2      63.53   -419.2
	layer straight   2    2    490.8739    -63.53     419.2      63.53    419.2
}

element TimoFourNodePlaneSec2D  1   1   2   3   4  5  1 
element TimoFourNodePlaneSec2D  2   4   5   6   7  5  1 
element TimoFourNodePlaneSec2D  3   7   8   9  10  5  1 
element TimoFourNodePlaneSec2D  4  10  11  12  13  5  1 
element TimoFourNodePlaneSec2D  5  13  14  15  16  5  1 
element TimoFourNodePlaneSec2D  7  16  17  18  19  5  1 

uniaxialMaterial  Elastic  3   1.e14 

node  101    -1.0   1219.2
fix 101 1 1 1

element truss  101  101  19  1.0  3

# ------------------------------------------------------

loadConst -time 0 

recorder Node -file TimoShear/react.out    -precision 12 -time -node  1  -dof 1 2 3 reaction;
recorder Node -file TimoShear/displ.out    -precision 12 -time -node  19  -dof 1 2 3 disp;


pattern Plain 2 Linear { 
    load  19  1.0e14   0.0   0.0
}

system BandGeneral
test NormDispIncr 1.0e-5 10
constraints Transformation
integrator LoadControl 0.01
algorithm Newton 
numberer RCM
analysis Static
set ok [analyze 1]


