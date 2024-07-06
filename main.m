clear
close all
% pi = sym(pi);
phi_semi = 60/180*pi;
xi_fov = 85/180*pi;
f = 1.5;
T = 1.0;
pt = 20;
rho_mi = 0.8;
ys = 2;
xd = 0; yd = 2; hd = 3;
% xs = 0.125; zs = 1.0212; 
xs = 0; zs = 0;
wm = 2.5/100; hm = 1.5/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
k = 30; l = 1;
coor_led = [0,ys,0].';
coor_PD = [xd,yd,hd].';
coor_RIS = [xs+wm/2+(l-1)*(wm+delta_wm), 0, zs+hm/2+(k-1)*(hm+delta_hm)].';

angle_led = [0,0];
angle_RIS = [pi/2,pi/2];
angle_PD = [pi,0];

p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);
% p_LoS = double(p_LoS);
coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
if length(coor_reflect) ~= 3
    p_NLoS = 0;
else
    p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
end
    % p_NLoS = double(p_NLoS);



