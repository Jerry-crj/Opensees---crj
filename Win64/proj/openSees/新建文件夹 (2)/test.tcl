
wipe

model basic -ndm 2 -ndf 2

node   1     0.0    0.0
node   2     1.0    1.0
node   3     2.0    0.0
node   4     1.0    2.0

fix 1 1 1 
fix 3 1 1  
fix 4 1 1 

mass 1  1  1
mass 2  1  1
mass 3  1  1
mass 4  1  1

uniaxialMaterial  Elastic     1  1.0e10

element  truss   1   1   2   1.0   1 
element  truss   2   3   2   1.0   1 
element  truss   3   4   2   1.0   1

timeSeries  Linear  2  -factor  1.0e10

pattern Plain 3 2 { 
    load   2   0.0   -1.0
}

# recorder   Pattern  pattern.out   -load   3    


# recorder  Node -file test.out -time  -node 4 -dof 1 2   disp
# recorder Node -file reactx.out -precision 12 -node 1 3 -dof 1 reaction;
# recorder Node -file reacty.out -precision 12 -node 1 3 -dof 2 reaction;
# recorder Node -file displx.out -precision 12 -node 2 -dof 1 disp;
# recorder Node -file disply.out -precision 12 -node 2 -dof 2 disp;


system BandGeneral
test NormDispIncr 1.0e-5 6 
constraints Transformation
integrator LoadControl  1.0e-4
algorithm Newton
numberer RCM
analysis Static
analyze 100

# wipe