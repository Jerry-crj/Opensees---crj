wipe

# 创建模型构造器
model basic -ndm 3 -ndf 3

# 定义节点
# 底部节点（固定）
node 1 0.0 0.0 0.0
node 2 1.0 0.0 0.0
node 3 1.0 1.0 0.0
node 4 0.0 1.0 0.0

# 顶部节点
node 5 0.0 0.0 1.0
node 6 1.0 0.0 1.0
node 7 1.0 1.0 1.0
node 8 0.0 1.0 1.0

# 顶部节点（加载位移）
node 12 1.0 0.0 0.0
node 13 1.0 1.0 0.0
node 16 1.0 0.0 1.0
node 17 1.0 1.0 1.0

# 定义边界条件
fix 1 1 1 1
fix 4 1 1 1
fix 5 1 1 1
fix 8 1 1 1
fix 12 1 1 1
fix 13 1 1 1
fix 16 1 1 1
fix 17 1 1 1

# 定义材料属性（假设弹性材料，杨氏模量=210000MPa，泊松比=0.3，密度=7850kg/m3）
uniaxialMaterial Elastic 10 1e12

uniaxialMaterial Concrete02 1 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 2 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 3 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 4 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 5 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 6 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 7 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 8 -40 -0.004 -5 -0.02 0.3 5 500
uniaxialMaterial Concrete02 9 -40 -0.004 -5 -0.02 0.3 5 500

nDMaterial MCTM 100 1 2 3 4 5 6 7 8 9
nDMaterial ElasticIsotropic 101 210000 0.3 7850

# 定义单元
element zeroLength 2 2 12 -mat 10 10 10 -dir 1 2 3
element zeroLength 3 3 13 -mat 10 10 10 -dir 1 2 3
element zeroLength 6 6 16 -mat 10 10 10 -dir 1 2 3
element zeroLength 7 7 17 -mat 10 10 10 -dir 1 2 3

element stdBrick 1 1 2 3 4 5 6 7 8 100

# 定义加载
timeSeries Linear 1 -factor 1.0e12

pattern Plain 1 1 {
    load 2 0.0 1.0 0.0
    load 3 0.0 1.0 0.0
    load 6 0.0 1.0 0.0
    load 7 0.0 1.0 0.0
}

recorder Node -file node1React.out -time -node 1 -dof 1 2 3 reaction
recorder Node -file node4React.out -time -node 4 -dof 1 2 3 reaction
recorder Node -file node5React.out -time -node 5 -dof 1 2 3 reaction
recorder Node -file node8React.out -time -node 8 -dof 1 2 3 reaction
recorder Node -file node2Displ.out -time -node 2 -dof 1 2 3 disp
recorder Node -file node3Displ.out -time -node 3 -dof 1 2 3 disp
recorder Node -file node6Displ.out -time -node 6 -dof 1 2 3 disp
recorder Node -file node7Displ.out -time -node 7 -dof 1 2 3 disp

recorder Element -file stress1.out -ele 1 material 1 stress
recorder Element -file stress2.out -ele 1 material 2 stress
recorder Element -file stress3.out -ele 1 material 3 stress
recorder Element -file stress4.out -ele 1 material 4 stress
recorder Element -file stress5.out -ele 1 material 5 stress
recorder Element -file stress6.out -ele 1 material 6 stress
recorder Element -file stress7.out -ele 1 material 7 stress
recorder Element -file stress8.out -ele 1 material 8 stress

recorder Element -file strain1.out -ele 1 material 1 strain
recorder Element -file strain2.out -ele 1 material 2 strain
recorder Element -file strain3.out -ele 1 material 3 strain
recorder Element -file strain4.out -ele 1 material 4 strain
recorder Element -file strain5.out -ele 1 material 5 strain
recorder Element -file strain6.out -ele 1 material 6 strain
recorder Element -file strain7.out -ele 1 material 7 strain
recorder Element -file strain8.out -ele 1 material 8 strain

recorder Element -file tstress1.out -ele 1 material 1 trussStress
recorder Element -file tstress2.out -ele 1 material 2 trussStress
recorder Element -file tstress3.out -ele 1 material 3 trussStress
recorder Element -file tstress4.out -ele 1 material 4 trussStress
recorder Element -file tstress5.out -ele 1 material 5 trussStress
recorder Element -file tstress6.out -ele 1 material 6 trussStress
recorder Element -file tstress7.out -ele 1 material 7 trussStress
recorder Element -file tstress8.out -ele 1 material 8 trussStress

recorder Element -file tstrain1.out -ele 1 material 1 trussStrain
recorder Element -file tstrain2.out -ele 1 material 2 trussStrain
recorder Element -file tstrain3.out -ele 1 material 3 trussStrain
recorder Element -file tstrain4.out -ele 1 material 4 trussStrain
recorder Element -file tstrain5.out -ele 1 material 5 trussStrain
recorder Element -file tstrain6.out -ele 1 material 6 trussStrain
recorder Element -file tstrain7.out -ele 1 material 7 trussStrain
recorder Element -file tstrain8.out -ele 1 material 8 trussStrain

# 定义分析参数
constraints Plain
numberer RCM
system BandGeneral
algorithm Newton
integrator LoadControl -1.0e-5 5
test NormDispIncr 1.0e-6 6
analysis Static

# 执行分析
analyze 2500


wipe