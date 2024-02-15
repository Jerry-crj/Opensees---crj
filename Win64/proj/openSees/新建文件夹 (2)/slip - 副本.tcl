
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
node 12 1.5 1.25
node 13 2.0 1.5
node 14 2.5 1.75
node 15 3.0 2.0
element PDSlip 11 1 2 3 4 11 2 2 2.0 1.0
element PDSlip 12 1 2 3 4 12 2 2 2.0 1.0
element PDSlip 13 1 2 3 4 13 2 2 2.0 1.0
element PDSlip 14 1 2 3 4 14 2 2 2.0 1.0
element PDSlip 15 1 2 3 4 15 2 2 2.0 1.0
element truss 16 11 12 1.0 3
element truss 17 12 13 1.0 3
element truss 18 13 14 1.0 3
element truss 19 14 15 1.0 3

node 21 1.0 2.0
node 22 1.5 1.75
node 23 2.0 1.5
node 24 2.5 1.25
node 25 3.0 1.0
element PDSlip 21 1 2 3 4 21 2 2 2.0 1.0
element PDSlip 22 1 2 3 4 22 2 2 2.0 1.0
element PDSlip 23 1 2 3 4 23 2 2 2.0 1.0
element PDSlip 24 1 2 3 4 24 2 2 2.0 1.0
element PDSlip 25 1 2 3 4 25 2 2 2.0 1.0
element truss 26 21 22 1.0 3
element truss 27 22 23 1.0 3
element truss 28 23 24 1.0 3
element truss 29 24 25 1.0 3


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
	load 4 1.0e4 2.0e4 
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