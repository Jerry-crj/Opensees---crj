
wipe 

model basic -ndm 2 -ndf 2

set AMatrix 1.0
set AFiber [expr 1.257e-3*5]

# material 1: fiber
uniaxialMaterial Hysteretic 11 1.0 0.05 0.001 1.80 0.1 40.0 -1.0 -0.05 -0.001 -0.80 -0.1 -40.0 0.8 0.2 0.0 0.0 
uniaxialMaterial Hysteretic 12 2.0 0.70 0.001 1.80 0.1 40.0 -2.0 -0.70 -0.001 -0.80 -0.1 -40.0 0.8 0.2 0.0 0.0 
uniaxialMaterial Parallel 1 11 12  
# material 2: rigid
uniaxialMaterial Elastic 2 1.0e10

# material 3: 断裂
uniaxialMaterial Hysteretic 3 1.0 0.05 0.001 0.10 0.1 40.0 -1.0 -0.05 -0.001 -0.80 -0.1 -40.0 0.8 0.2 0.0 0.0 



# uniaxialMaterial Concrete02 2 -40 -0.0025 -8 -0.015 0.3 1.4 1000
# uniaxialMaterial Hysteretic 3 1500.0 0.04 1.0 0.08 10.0 10.0 -1.0 -1.0 -2.0 -2.0 -3.0 -3.0 0.8 0.2 0.0 0.0 

nDMaterial ElasticIsotropic 10 1.0e7 0.25

node 1 0.0 0.0
node 2 1.0 0.0
node 3 0.0 1.0
node 4 1.0 1.0
node 5 0.0 1.0
node 6 1.0 1.0
node 7 0.0 2.0
node 8 1.0 2.0
node 107 0.0 2.0
node 108 1.0 2.0
node 1001 0.5 0.5
node 1002 0.5 1.5

fix 1 1 1
fix 2 1 1
fix 107 1 1
fix 108 1 1

element quad 1 1 2 4 3 1.0 "PlaneStrain" 10
element quad 2 5 6 8 7 1.0 "PlaneStrain" 10
element truss 3 1001 1002 1.0 1
element PDSlip2D 11 1 2 4 3 1001 2 2 0.0 1.0
element PDSlip2D 12 5 6 8 7 1002 2 2 0.0 1.0
element zeroLength 103 3 5 -mat 2 3 -dir 1 2
element zeroLength 104 4 6 -mat 2 3 -dir 1 2
element zeroLength 107 107 7 -mat 2 2 -dir 1 2
element zeroLength 108 108 8 -mat 2 2 -dir 1 2

recorder Node -file displ.out -precision 8 -node 7 8 -dof 2 disp
recorder Node -file react.out -precision 8 -node 1 2 -dof 2 reaction
recorder Element -file ele.out -precision 8 -ele 3 material stressStrain


timeSeries Linear 1 -factor 1.0e10
pattern Plain 1 1 {
	load 7 0.0 1.0
	load 8 0.0 1.0
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
analyze 1500


wipe


