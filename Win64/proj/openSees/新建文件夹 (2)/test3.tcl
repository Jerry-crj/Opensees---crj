# 纯压

wipe 
model basic -ndm 2 -ndf 2 

node  1   0.0   0.0 
node  2   1.0   0.0 
node  3   0.0   1.0 
node  4   1.0   1.0 
node  5   1.0   1.0 
node  6   1.0   1.0 

fix  1  1  1 
fix  2  0  1 

equalDOF     5  4  1  2 
equalDOF     6  5  1  2 



uniaxialMaterial  Elastic  1   1.0e5


element truss  1  1  2  1.0   1
element truss  2  1  3  1.0   1
element truss  3  1  4  1.0   1
element truss  4  2  3  1.0   1
element truss  5  2  4  1.0   1
element truss  6  3  4  1.0   1

recorder Node -file nodedisp1.out -precision 16   -node  1  -dof 1 2 disp
recorder Node -file nodedisp2.out -precision 16   -node  2  -dof 1 2 disp
recorder Node -file nodedisp3.out -precision 16   -node  3  -dof 1 2 disp
recorder Node -file nodedisp4.out -precision 16   -node  4  -dof 1 2 disp
recorder Node -file nodedisp5.out -precision 16   -node  5  -dof 1 2 disp

pattern Plain 1 Linear {
	load   3   0.0   1.0e3    
	load   5   0.0   1.0e3    
}

constraints Transformation;     				
numberer RCM;					
system BandGeneral; 			
test NormDispIncr  0.01  10  1			
algorithm Newton;					
integrator LoadControl  0.1
analysis Static	
analyze 1

wipe