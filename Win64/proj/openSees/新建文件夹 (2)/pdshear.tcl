# 纯压

wipe 
model basic -ndm 2 -ndf 2 

node  1      0.0      0.0 
node  2   1000.0   1000.0 

fix  1  1  1 

uniaxialMaterial Elastic  1   200
uniaxialMaterial Elastic  2   100
nDMaterial  CombinedUniaxialTimoMat  3   1  2

element PDTrussWithShear2D  1  1  2  1.0   3 
# element truss  1  1  2  1.0   1

recorder Node -file react.out  -node 1  2 -dof 1 2 reaction;
recorder Node -file displ.out  -node 1  2 -dof 1 2 disp;

timeSeries Linear 1

pattern Plain 1 1 { 
   load  2  1.0  1.0 
}

system BandGeneral
test NormDispIncr 1.e-8 6
constraints Transformation
integrator LoadControl  1.0
algorithm Newton 
numberer RCM
analysis Static
analyze 1

wipe


