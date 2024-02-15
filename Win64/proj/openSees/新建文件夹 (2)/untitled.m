clear; close all; clc

x0 = 2; 
y0 = 1;
a = 0.5;
b = 0.5;
c = -1;

x = -(a * c - b^2 * x0 + a* b* y0)/(a^2 + b^2);
y = -( b* c + a* b* x0 - a^2 * y0)/(a^2 + b^2);

