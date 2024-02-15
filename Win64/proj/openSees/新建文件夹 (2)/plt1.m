clear; close all; clc

load("stress.out")
load("strain.out")

figure
plot(strain(:,2),stress(:,2))
title("正应力-正应变")

figure
plot(strain(:,3),stress(:,3))
title("剪应力-剪应变")
