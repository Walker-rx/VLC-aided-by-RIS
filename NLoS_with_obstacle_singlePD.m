clear
close all

phi_semi = 60/180*pi;
xi_fov = 85/180*pi;
f = 1.5; 
T = 1.0;
pt = 20;
rho_mi = 0.8;

B = 200e6;
q = 3;
R_pd = 0.53;
N = 1e-13;

xs=3; ys=3; zs=0;
led_polar = 0; led_azi = 0;
coor_led = [xs,ys,zs].';
angle_led = [led_polar,led_azi];

wm = 4/100; hm = 4/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
total_length = 6; total_hight = 1;
xr_max = 6:0.2:8; zr_max = 1:0.1:3;
max_k = total_hight/hm;    max_l=total_length/wm;

k = 8; l = 20;
delta_wm = round( (((100*total_length)-l*100*wm)/(l-1))*10 )/10 /100;
delta_hm = round( (((100*total_hight)-k*100*wm)/(k-1))*10 )/10 /100;
if delta_wm <1/100
    delta_wm = 1/100;
    l = floor((total_length-wm)/(wm+delta_wm))+1;
end
if delta_hm <1/100
    delta_hm = 1/100;
    k = floor((total_hight-hm)/(hm+delta_hm))+1;
end

% obstacle_length = 2.5;
obstacle_length = 2.5;
coor_obstacle = [3,3,1];
% angle_obstacle = [pi 0];
angle_obstacle = [pi/4 pi/4];

xd_save = 1:0.5:7;
yd_save = 1:0.5:7;
% xd_for_plot = 9;  zd_for_plot = 9;
xd_for_count = 6;  yd_for_count = 6;
hd = 3;
pd_polar = pi/3+pi/2; pd_zai = pi/3;
% pd_polar = pi; pd_zai =0;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
coor_PD = [xd_save(xd_for_count),yd_save(yd_for_count),hd].';
angle_PD = [pd_polar,pd_zai];

p_NLoS_save = [];
rate_LoS_save = [];
rate_total_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

if ~judge_intersect(coor_led,coor_PD,coor_obstacle,angle_obstacle,obstacle_length)
    p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T); 
else
    p_LoS = 0;
end
rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);

for i = 1:length(xr_max)
    for j = 1:length(zr_max)
        xr = xr_max(i)-total_length;   zr = zr_max(j)-total_hight;
        % str_diff=strcat('diff_'+string(i)+'_'+string(j));
        % eval(strcat(str_diff+'=[];'));
        str_x_use=strcat('RIS_x_use_'+string(i)+'_'+string(j));
        eval(strcat(str_x_use+'=[];'));
        str_z_use=strcat('RIS_z_use_'+string(i)+'_'+string(j));
        eval(strcat(str_z_use+'=[];'));
        str_x_unuse=strcat('RIS_x_unuse_'+string(i)+'_'+string(j));
        eval(strcat(str_x_unuse+'=[];'));
        str_z_unuse=strcat('RIS_z_unuse_'+string(i)+'_'+string(j));
        eval(strcat(str_z_unuse+'=[];'));

        p_NLoS_singlePD = 0;
        
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);
                if all((angle_RIS > 0) & (angle_RIS < pi))
                    coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                    if all(abs(coor_reflect-coor_RIS)<1e-7)
                        if (~judge_intersect(coor_led,coor_reflect,coor_obstacle,angle_obstacle,obstacle_length)) && (~judge_intersect(coor_reflect,coor_PD,coor_obstacle,angle_obstacle,obstacle_length))
                            p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        else
                            p_NLoS = 0;
                            fprintf('NLoS obstacle happens in i=%d ,j=%d ,ii=%d ,jj=%d \n',i,j,ii,jj);
                        end
                                         
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;
                        if p_NLoS == 0
                            eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                        else
                            eval(strcat(str_x_use,'=[',str_x_use,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_use,'=[',str_z_use,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));                            
                        end                       
                    else
                        diff_tmp = NaN;
                        fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
                        eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                        eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                        eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                    end                    
                else
                    diff_tmp = NaN;
                    fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
                    eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)));
                    eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                    eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                end
            end
        end

        rate_total = minrate(B,pt,q,R_pd,N,(p_LoS+p_NLoS_singlePD)/pt);
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
        rate_total_save = [rate_total_save rate_total];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(unique(xr_max-total_length/2), unique(zr_max-total_hight/2));
Z_NLoS_power = reshape(p_NLoS_save,length(zr_max),length(xr_max));
figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_NLoS_power(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('NLOS功率');
% title('NLOS功率');
% 
% figure
% contour(X, Y, Z_NLoS_power, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Z');
% title('NLOS功率');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_rate = reshape(rate_total_save,length(zr_max),length(xr_max));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_rate(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Z');
% zlabel('NLOS速率');
% title('NLOS速率');

figure
[C,h]=contour(X, Y, Z_total_rate, 6);
clabel(C,h)
xlabel('X');
ylabel('Z');
title('NLOS速率');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% globalMin = min(min(Z_total_power(:)), min(Z_LoS_power(:)));
% globalMax = max(max(Z_total_power(:)), max(Z_LoS_power(:)));
% 
% subplot(2, 1, 1);
% [C,h]=contour(X, Y, Z_total_power);
% clabel(C,h)
% xlabel('X');
% ylabel('Y');
% title('功率总和');
% 
% subplot(2, 1, 2);
% [C,h]=contour(X, Y, Z_LoS_power,6);
% clabel(C,h)
% xlabel('X');
% ylabel('Y');
% title('LOS功率');
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% globalMin = min(min(Z_total_rate(:)), min(Z_LoS_rate(:)));
% globalMax = max(max(Z_total_rate(:)), max(Z_LoS_rate(:)));
% 
% subplot(2, 1, 1);
% [C,h]=contour(X, Y, Z_total_rate);
% clabel(C,h)
% xlabel('X');
% ylabel('Y');
% title('存在NLoS径时的最低速率');
% 
% subplot(2, 1, 2);
% [C,h]=contour(X, Y, Z_LoS_rate,6);
% clabel(C,h)
% xlabel('X');
% ylabel('Y');
% title('只有LoS径时的最低速率');

