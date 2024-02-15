
wipe 

model basic -ndm 2 -ndf 2

uniaxialMaterial Hysteretic 11 1.0 0.05 0.001 1.80 0.1 40.0 -1.0 -0.05 -0.001 -0.80 -0.1 -40.0 0.8 0.2 0.0 0.0 
uniaxialMaterial Hysteretic 12 2.0 0.70 0.001 1.80 0.1 40.0 -2.0 -0.70 -0.001 -0.80 -0.1 -40.0 0.8 0.2 0.0 0.0 
uniaxialMaterial Parallel 1 11 12 
uniaxialMaterial Elastic 2 1.0e10

node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 0.0

fix 1 1 1
fix 2 0 1
fix 3 1 1

element truss 1 1 2 1.0 1
element zeroLength 2 2 3 -mat 2 -dir 1

recorder Node -file displ.out -precision 8 -time -node 2 -dof 1 disp
recorder Node -file react.out -precision 8 -time -node 1 -dof 1 reaction

timeSeries Linear 1 -factor 1.0e10
pattern Plain 1 1 {
	load 2 1.0 0.0
}

set dt  1.0e-4
set time [getTime]
set tol 3.0e-4

constraints Transformation;     				
numberer Plain
system BandGeneral		
test NormDispIncr  $tol  10	
algorithm Newton;					
integrator LoadControl  $dt
analysis Static	
analyze 20000

wipe


