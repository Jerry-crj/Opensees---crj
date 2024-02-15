# 纯压

wipe 

model basic -ndm 2 -ndf 3 

node  1      0.0      0.0 
node  2   1800.0      0.0 

fix  1  1  1  1

uniaxialMaterial  Elastic  1  2.0e5

section Fiber 1  {
	patch  rect  1   20   1  -100  -50   100   50
}

geomTransf  Linear  1

element dispBeamColumn  1  1  2  5  1  1

recorder Node -file ereact.out -time -node 1 -dof 1 2 3 reaction;
recorder Node -file edispl.out -time -node 2 -dof 1 2 3 disp;

timeSeries Linear 1

pattern Plain 1 1 { 
   load  2  0.0  1.0  0.0
}

system BandGeneral
test NormDispIncr 1.e-8 6
constraints Transformation
integrator DisplacementControl  2  2  0.005
algorithm Newton 
numberer RCM
analysis Static
analyze 10


wipe


