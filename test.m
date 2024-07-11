% clear 
% close all


thresholds = [0 0.25 4];
figure;

custommap = [
    0 0 1
    0 1 0
    1 0 0];
[C, ~] = contourf(X, Y, Z_diff, thresholds);
colormap(custommap)

clabel(C)
