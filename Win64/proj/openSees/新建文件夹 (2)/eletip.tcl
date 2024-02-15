# 纯压

wipe 

model basic -ndm 2 -ndf 3 

node  1      0.0      0.0 
node  2    600.0      0.0 
node  3   1200.0      0.0 
node  4   1800.0      0.0 

fix  1  1  1  1

uniaxialMaterial  Elastic  1  2.0e5
uniaxialMaterial  Elastic  2  8.0e10
nDMaterial  CombinedUniaxialTimoMat  3   1  2

section FiberSectionTimoConstShear 1  {
	patch  rect  3   20   1  -100  -50   100   50
}

element TimoFourNodePlaneSec2D  1  1  2  3  4  5  1

recorder Node -file treact.out -time -node 1 -dof 1 2 3 reaction;
recorder Node -file tdispl.out -time -node 4 -dof 1 2 3 disp;

timeSeries Linear 1

pattern Plain 1 1 { 
   load  4  0.0  1.0  0.0
}

system BandGeneral
test NormDispIncr 1.e-8 6 
constraints Transformation
integrator DisplacementControl  4  2  0.005
algorithm Newton 
numberer RCM
analysis Static
analyze 10


wipe


