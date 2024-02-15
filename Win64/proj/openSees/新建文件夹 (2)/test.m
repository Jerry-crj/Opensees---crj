clear; close all; clc

load("TimoShear\displ.out");
load("TimoShear\react.out");

figure
plot(displ(:,2),-react(:,2),'r');


