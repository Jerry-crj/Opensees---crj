clear; close all; clc

load node2Displ.out;
load node3Displ.out;
load node6Displ.out;
load node7Displ.out;

load node1React.out;
load node4React.out;
load node5React.out;
load node8React.out;

displ = node2Displ+node3Displ+node6Displ+node7Displ;
react = node1React+node4React+node5React+node8React;
displ = displ/4;
react = -react;

figure 
plot(displ(:,3),react(:,3));


