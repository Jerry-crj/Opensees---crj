clear; close all; clc

load displ.out
load react.out

figure
plot(sum(displ,2),-sum(react,2))

load ele.out
figure
plot(ele(:,2),ele(:,1))