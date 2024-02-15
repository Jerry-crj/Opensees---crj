clear; close all; clc

load("stress.out");
load("strain.out");
load("dmg.out")

figure
plot(strain(:,1),stress(:,1))
grid on

figure
plot(strain(:,2),stress(:,2))
grid on