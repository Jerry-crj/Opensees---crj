# 纯压

wipe 
model basic -ndm 2 -ndf 2 

node  1      0.0      0.0 
node  2   1000.0      0.0 
node  3      0.0   1000.0 
node  4   1000.0   1000.0 

mass  3   1000   1000
mass  4   1000   1000

fix  1  1  1 
fix  2  1  1 

equalDOF 3 4 1 2

nDMaterial ovalConcrete  1  -40   -2.e-3   -10.  -1.e-2    0.1    2.   2e4   2.   0.4  8000   
# nDMaterial ovalConcrete  3  -14.0   -0.0014    -2.5   -0.01    0.1   2.8    10000   2.8    0.4   4000

nDMaterial TimoToQuadMat  2  1    

element quad  1  1  2  4  3  1.0  "PlaneStrain"  2 

uniaxialMaterial  Elastic  3   1.0e10

node 11     0.0   1100.0
node 12  1000.0   1100.0
node 13  -100.0   1000.0
node 14  1100.0   1000.0

fix 11 1 1
fix 12 1 1
fix 13 1 1
fix 14 1 1

element  truss   11   11  3   1.0   3
element  truss   12   12  4   1.0   3
element  truss   13   13  3   1.0   3
element  truss   14   14  4   1.0   3

puts "model has been built..." 


recorder Node -file dd.out -time -node 3 4 -dof 1 2 disp;
recorder Element -file stress.out -time -ele 1 material 1 stress
recorder Element -file strain.out -time -ele 1 material 1 strain
recorder Element -file dam.out -time -ele 1 material 1 damage

timeSeries Linear 1

pattern Plain 1 1 { 
   load  4  1.0  1.0 
   load  3  1.0  1.0 
}

system BandGeneral
test NormDispIncr 1.e-8 6
constraints Transformation
integrator LoadControl 1.0e5
algorithm Newton 
numberer RCM
analysis Static
analyze  200

# loadConst  -time 0
# pattern Plain 2 1 { 
   # load  4  1.0  0.0 
   # load  3  1.0  0.0 
# }

test NormDispIncr 1.e-5 6  
integrator LoadControl -1.e5
analyze 400

test NormDispIncr 1.e-5 6  
integrator LoadControl 1.e5
analyze 600

test NormDispIncr 1.e-5 6  
integrator LoadControl -1.e5
analyze 800



# test NormDispIncr 1.e-5 6  
# integrator LoadControl -5.e4
# analyze 400


# loadConst  -time 0
# pattern Plain 3 1 { 
   # load  4  0.0  1.0 
   # load  3  0.0  1.0 
# }
# integrator LoadControl 5.e5
# analyze 1510



wipe


