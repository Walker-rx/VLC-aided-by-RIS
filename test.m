clear 
close all

% 创建网格化数据
[X, Y] = meshgrid(1:5, 1:5);
Z = X .^ 2 + Y .^ 2;

% 绘制等高线图
figure;
contourf(X, Y, Z,ShowText="on");

% 添加数值标注
clabel(Z);

% 自定义图形
title('Contour Plot with Labels');
xlabel('X-axis');
ylabel('Y-axis');