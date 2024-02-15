clear; close all; clc

load displ.out
load react.out

figure
plot(displ(:,2),-react(:,2))
