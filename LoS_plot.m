clear
close all

phi_semi = 60/180*pi;
xi_fov = 85/180*pi;
f = 1.5;
T = 1.0;
pt = 20;
rho_mi = 0.8;

xs=3; ys=3; zs=0;
led_polar = 0; led_azi = 0;
coor_led = [xs,ys,zs].';
angle_led = [led_polar,led_azi];

hd = 3;
pd_polar = pi; pd_zai = 0;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
angle_PD = [pd_polar,pd_zai];

xd_save = 1:0.5:7;
yd_save = 1:0.5:7;
p_LoS_save = [];
for i = 1:length(xd_save)
    xd = xd_save(i);
    for j = 1:length(yd_save)
        yd = yd_save(j);
        coor_PD = [xd,yd,hd].'; 
        p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);
        p_LoS_save = [p_LoS_save p_LoS];                
    end
end
[X, Y] = meshgrid(unique(xd_save), unique(yd_save));
Z = reshape(p_LoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('直达径功率');








