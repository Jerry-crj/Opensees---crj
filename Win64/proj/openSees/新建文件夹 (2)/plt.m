clear; close all; clc

load("reactx.out");
load("reacty.out");
load("displx.out");
load("disply.out");

figure
plot(disply(:,1),-(reacty(:,1)+reacty(:,2)))

