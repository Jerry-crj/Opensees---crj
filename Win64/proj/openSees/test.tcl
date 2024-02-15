# 清空所有先前的模型定义
wipe

# 创建模型
model basic -ndm 2 -ndf 3

# 定义节点
node 1 0.0 0.0
node 2 0.0 1000.0
fix 1 1 1 1
mass 2 1.0 1.0 0.0

# 定义元素
geomTransf Linear 1
element elasticBeamColumn 1 1 2 1.0 1.0 1.0 1

# 定义推覆载荷

timeSeries Linear 1 -factor 1.0
pattern Plain 1 1 {
    load 2 1.0e5 0.0 0.0
}

# 设置分析选项
constraints Plain
numberer Plain
system FullGeneral
test NormDispIncr 1.0e-5 10
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient
# integrator LoadControl 0.1
# analysis Static

# 进行分析
for {set i 0} {$i < 10} {incr i} {
analyze 1 0.01
reactions
puts "[nodeReaction 2 1]"
}


wipe