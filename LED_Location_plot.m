clear
close all

phi_semi = 60/180*pi;
xi_fov = 85/180*pi;
f = 1.5; 
T = 1.0;
pt = 20;
rho_mi = 0.8;

xd = 5; yd = 1; hd = 3;
coor_RIS = [5, 0, hd/2].';
angle_RIS = [pi/2,pi/2];
coor_PD = [xd,yd,hd].';
angle_PD = [pi,0];

wm = 2.5/100; hm = 1.5/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
wd = 0.1/100; ld = 0.1/100;
A_PD = wd*ld;
angle_led = [0,0];

reflect_point_save = [];
p_NLoS_save = [];
p_LoS_save = [];
xs_save = [];
ys_save = [];
diff_save = [];

% xs = 0.125; zs = 1.0212; 
% xs = 0; zs = 0;
% k = 30; l = 1;
% coor_led = [0,ys,0].';
% xd_save = -4:0.5:4;
% yd_save = 0.5:0.5:8;
% coor_RIS = [xs+wm/2+(l-1)*(wm+delta_wm), 0, zs+hm/2+(k-1)*(hm+delta_hm)].';
% coor_RIS = [0, 0, hd/2].';
% angle_led = [0,0];
% angle_RIS = [pi/2,pi/2];
% reflect_point_save = [];
loop = 0;
for xs = 3.5:0.001:6.5
    for ys = 0.5:0.001:3
        coor_led = [xs,ys,0].';
        coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
        if isempty(coor_reflect)
            continue;
        else
            reflect_point_save = [reflect_point_save;coor_reflect];
            p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
            p_NLoS_save = [p_NLoS_save p_NLoS];
            xs_save = [xs_save,xs];
            ys_save = [ys_save,ys];
            loop = loop+1;

            p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);
            p_LoS_save = [p_LoS_save p_LoS];

            diff_save = [diff_save,p_LoS/p_NLoS];
        end
    end
end
loop
a = [xs_save;ys_save;zeros(1,length(xs_save))];
plot3(a(1,:), a(2,:), a(3,:), 'o');

[X, Y] = meshgrid(unique(xd_save), unique(yd_save));
% Z_LoS = reshape(p_LoS_save,length(yd_save),length(xd_save));
Z_NLoS = reshape(p_NLoS_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% title("LoS")
% xlabel('X');
% ylabel('Y');
% zlabel('LOS功率');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_NLoS(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
title("NLoS")
xlabel('极角');
ylabel('方位角');
zlabel('NLOS功率');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_distance = reshape(distance_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_distance(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
title("距离")
xlabel('极角');
ylabel('方位角');
zlabel('距离');

% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS(:)+Z_NLoS(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% title("LoS+NLoS")
% xlabel('X');
% ylabel('Y');
% zlabel('LoS+NLOS功率');

% figure;
% scatter3(reflect_point_save(:,1),reflect_point_save(:,2),reflect_point_save(:,3))
% title("反射点位置")
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% reflect_point_xmax = max(reflect_point_save(:,1));
% reflect_point_xmin = min(reflect_point_save(:,1));
% reflect_point_zmax = max(reflect_point_save(:,3));
% reflect_point_zmin = min(reflect_point_save(:,3));








