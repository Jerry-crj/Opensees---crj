clear; close all; clc

load node5Displ.out;
load node6Displ.out;
load node7Displ.out;
load node8Displ.out;

load node1React.out;
load node2React.out;
load node3React.out;
load node4React.out;

displ = node5Displ+node6Displ+node7Displ+node8Displ;
react = node1React+node2React+node3React+node4React;
displ = displ/4;
react = -react;

figure 
plot(displ(:,3),react(:,3));


