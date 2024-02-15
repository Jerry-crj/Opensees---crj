# 纯压

wipe 
model basic -ndm 2 -ndf 2 

node  1   0.0   0.0 
node  2   1.0   0.0 
node  3   0.0   1.0 
node  4   1.0   1.0 

fix  1  1  1 
fix  2  1  1 

equalDOF 3 4 1 2

nDMaterial  ovalConcrete  1  -20.0   -0.002    -4.   -0.01    0.05   2.0    500   2.0   0.1   200  
nDMaterial TimoToQuadMat  2  1    

element quad  1  1  2  4  3  1.0  "PlaneStrain"  2 

uniaxialMaterial  Elastic  3   1.0e14

node 11   0.0   2.0
node 12   1.0   2.0
node 13  -1.0   1.0
node 14   2.0   1.0

fix 11 1 1
fix 12 1 1
fix 13 1 1
fix 14 1 1

element  truss   11   11  3   1.0   3
element  truss   12   12  4   1.0   3
element  truss   13   13  3   1.0   3
element  truss   14   14  4   1.0   3

puts "model has been built..." 


recorder Node    -file displ.out   -node 3 4 -dof 1 2 disp;
recorder Element -file stress.out  -ele  1  material 1 stress
recorder Element -file strain.out  -ele  1  material 1 strain
recorder Element -file dmg.out     -ele  1  material 1 dmgfct
# recorder Element -file stiff.out   -ele  1  material 1 tangent

timeSeries Linear 1

pattern Plain 1 1 { 
   load  4  0.0   1.0
   load  3  0.0   1.0
}

system BandGeneral
test NormDispIncr 1.e-8 6
constraints Transformation
integrator LoadControl  1.e-5
algorithm Newton 
numberer RCM
analysis Static
analyze 1

loadConst -time 0

pattern Plain 2 1 { 
   load  4   1.0e14    0.0 
   load  3   1.0e14    0.0 
}

integrator LoadControl 1.e-5
analyze 1500

integrator LoadControl  1.e-5  1
integrator LoadControl -1.e-5
analyze 2000

# integrator LoadControl -1.e-5
# analyze 1200

wipe


