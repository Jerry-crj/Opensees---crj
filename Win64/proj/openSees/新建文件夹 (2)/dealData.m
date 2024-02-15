clear; close all; clc

load("TimoShear\displ.out");

rst = displ(:,2);
for i = 1:length(rst)
    rstround(i) = round(rst(i),3);
end

fid = fopen("input.in","w");
for i = 1:length(rst)
    fprintf(fid,"%11.4f \n", rstround(i));
end
fclose("all");