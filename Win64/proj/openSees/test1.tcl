wipe
# 定义模型参数
model basic -ndm 2 -ndf 2

set dt 1.0e-5

# 创建节点
node 1 0.0 0.0
node 2 1.0 0.0
node 3 1.0 0.0

# 创建节点质量
mass 1 1.0 1.0
mass 2 1.0 1.0
mass 3 1.0 1.0

# 设置节点约束
fix 1 1 1
fix 2 0 1
fix 3 1 1

# 定义材料属性
uniaxialMaterial Elastic 10 1.0e10
uniaxialMaterial Concrete02 1 -32.594e3 -0.003243 -6.5188e3 -0.0489 0.75 3.2594e3 1100

# 定义弹簧杆单元
element truss 1 1 2 1.0 1
element zeroLength 2 2 3 -mat 10 10 -dir 1 2 

# 应用荷载
timeSeries Linear 1 -factor 1.0e10
pattern Plain 1 1 {
    load 2 1.0 0.0
}

recorder Node -file displ.out -time -node 2 -dof 1 disp;
recorder Node -file react.out -time -node 1 -dof 1 reaction;

# 定义分析参数
# constraints Penalty 1.0e8 1.0e0
constraints Plain
numberer Plain
system FullGeneral
test NormUnbalance 1.0e-6 10 0
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient

# 执行分析
for {set i 1} {$i <= 10} {incr i} {
    analyze 1 0.001
    puts "[nodeUnbalance 2 1] [nodeUnbalance 2 2]"
}

# # 执行分析
# for {set i 1} {$i <= 10} {incr i} {
    # remove timeSeries 1
    # remove loadPattern 1
    # timeSeries Linear 1 -factor [expr $i*1.0e10]
    # pattern Plain 1 1 {
        # load 2 1.0 0.0
    # }
    # if {$i <= 9} {
    # analyze 1
    # } else {
        # analyze 10
    # }
# }

wipe