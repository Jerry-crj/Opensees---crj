
wipe
# 定义模型
model basic -ndm 2 -ndf 2

# 定义材料属性
uniaxialMaterial Elastic 1 2.0 
uniaxialMaterial Elastic 2 1.0e10 
uniaxialMaterial Elastic 3 2.0e5 

# 定义节点
node 1 0.0 0.0
node 2 4.0 0.0
node 3 4.0 4.0
node 4 0.0 4.0

# 定义单元
element truss 1 1 2 1.0 3
element truss 2 1 3 1.0 3
element truss 3 1 4 1.0 3
element truss 4 2 3 1.0 3
element truss 5 2 4 1.0 3
element truss 6 3 4 1.0 3

node 11 1.0 1.0
element PDSlip 11 1 2 3 4 11 2 2 2.0 1.0

# 定义边界条件
fix 1 1 1
fix 2 1 1
# fix 3 1 1
# fix 4 1 1

# 定义recorder
recorder Node -file "displacement.out" -time -node 1 2 3 4 5 6 -dof 1 2 disp
recorder Element -file "ele.out" -time -ele 10 deform

# 定义荷载
pattern Plain 1 Linear {
	load 11 1.0e4 2.0e4 
#	0.0
}

# 执行分析
system BandGeneral
numberer RCM
constraints Plain
algorithm Newton
integrator LoadControl 1.0
test NormDispIncr 1.0e-6 5 1
analysis Static

analyze 3

# print A
wipe