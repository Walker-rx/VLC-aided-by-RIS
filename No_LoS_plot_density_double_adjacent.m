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

hd = 1.8;
pd_polar = pi/4; pd_zai = -pi/4;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
angle_PD = [pd_polar,pd_zai];

wm = 4/100; hm = 4/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
xr = 1;   zr = 2; yr = 1;
total_length = 7-xr;  total_hight = 3-zr;
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
RIS_x = xr+wm/2:(wm+delta_wm):xr+wm/2+(l-1)*(wm+delta_wm);
RIS_z = zr+hm/2:(hm+delta_hm):zr+hm/2+(k-1)*(hm+delta_hm);

xd_save = 1:0.5:7;
yd_save = 1:0.5:7;
xd_for_plot = 9;  yd_for_plot = 9;
xd_for_plot_yz = 9;  yd_for_plot_yz = 9;

p_LoS_save = [];
p_NLoS_save_xz = [];
p_NLoS_save_yz = [];
p_NLoS_save = [];

rate_LoS_save = [];
rate_NLoS_save_xz = [];
rate_NLoS_save_yz = [];
rate_total_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

% coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(xd_save)
    xd = xd_save(i);
    for j = 1:length(yd_save)
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
        
        % str_diff_yz=strcat('diff_yz_'+string(i)+'_'+string(j));
        % eval(strcat(str_diff_yz+'=[];'));
        str_y_use_yz=strcat('RIS_y_use_yz_'+string(i)+'_'+string(j));
        eval(strcat(str_y_use_yz+'=[];'));
        str_z_use_yz=strcat('RIS_z_use_yz_'+string(i)+'_'+string(j));
        eval(strcat(str_z_use_yz+'=[];'));
        str_y_unuse_yz=strcat('RIS_y_unuse_yz_'+string(i)+'_'+string(j));
        eval(strcat(str_y_unuse_yz+'=[];'));
        str_z_unuse_yz=strcat('RIS_z_unuse_yz_'+string(i)+'_'+string(j));
        eval(strcat(str_z_unuse_yz+'=[];'));
        
        yd = yd_save(j);
        coor_PD = [xd,yd,hd].';
        p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);  
        % if p_LoS ~= 0
        %     fprintf("LoS exists in i=%d,j=%d \n",i,j);
        % end
        p_NLoS_singlePD = 0;
        p_NLoS_xz_tmp = 0;
        p_NLoS_yz_tmp = 0;
        rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);

                coor_RIS_yz = [0,yr+wm/2+(jj-1)*(wm+delta_wm), zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS_yz = findReflect_angle_yz(coor_led,coor_PD,coor_RIS_yz);
                if all((angle_RIS > 0) & (angle_RIS < pi))
                    coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                    if all(abs(coor_reflect-coor_RIS)<1e-7)
                        p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        % diff_tmp = p_LoS/p_NLoS;
                        % eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;
                        p_NLoS_xz_tmp = p_NLoS_xz_tmp + p_NLoS;
                        if p_NLoS == 0
                            eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                        else
                            eval(strcat(str_x_use,'=[',str_x_use,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_use,'=[',str_z_use,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));                            
                        end                       
                    else
                        % diff_tmp = NaN;
                        fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
                        % eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                        eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                        eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                    end                    
                else
                    % diff_tmp = NaN;
                    fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
                    % eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)));
                    eval(strcat(str_x_unuse,'=[',str_x_unuse,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                    eval(strcat(str_z_unuse,'=[',str_z_unuse,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% yz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (angle_RIS_yz(1) > 0) && (angle_RIS_yz(1) < pi) && (angle_RIS_yz(2) > -pi/2) && (angle_RIS_yz(2) < pi/2)
                    coor_reflect_yz = findReflect_point(coor_led,coor_PD,coor_RIS_yz,angle_RIS_yz,wm,hm);
                    if all(abs(coor_reflect_yz-coor_RIS_yz)<1e-7)
                        p_NLoS_yz = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect_yz,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        % diff_tmp_yz = p_LoS/p_NLoS_yz;
                        % eval(strcat(str_diff_yz+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp_yz)+';'));
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS_yz;
                        p_NLoS_yz_tmp = p_NLoS_yz_tmp + p_NLoS_yz;
                        if p_NLoS_yz == 0
                            eval(strcat(str_y_unuse_yz,'=[',str_y_unuse_yz,','+string(yr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_unuse_yz,'=[',str_z_unuse_yz,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                        else
                            eval(strcat(str_y_use_yz,'=[',str_y_use_yz,','+string(yr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_use_yz,'=[',str_z_use_yz,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));                            
                        end                       
                    else
                        % diff_tmp_yz = NaN;
                        fprintf("NaN happen in yz wall, i=%d,j=%d,ii=%d,jj=%d \n",i,j,ii,jj);
                        % eval(strcat(str_diff_yz+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp_yz)+';'));
                        eval(strcat(str_y_unuse_yz,'=[',str_y_unuse_yz,','+string(yr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                        eval(strcat(str_z_unuse_yz,'=[',str_z_unuse_yz,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                    end                    
                else
                    % diff_tmp_yz = NaN;
                    fprintf("NaN happen in yz wall, i=%d,j=%d,ii=%d,jj=%d \n",i,j,ii,jj);
                    % eval(strcat(str_diff_yz+'('+string(ii)+','+string(jj)+')'+'=NaN;'));
                    eval(strcat(str_y_unuse_yz,'=[',str_y_unuse_yz,','+string(yr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                    eval(strcat(str_z_unuse_yz,'=[',str_z_unuse_yz,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                end

                
            end
        end

        rate_total = minrate(B,pt,q,R_pd,N,(p_LoS+p_NLoS_singlePD)/pt);
        p_LoS_save = [p_LoS_save p_LoS];
        p_NLoS_save_xz = [p_NLoS_save_xz p_NLoS_xz_tmp];
        p_NLoS_save_yz = [p_NLoS_save_yz p_NLoS_yz_tmp];
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
                             
        rate_LoS_save = [rate_LoS_save rate_LoS];
        rate_NLoS_save_xz = [rate_NLoS_save_xz minrate(B,pt,q,R_pd,N,p_NLoS_xz_tmp/pt)];
        rate_NLoS_save_yz = [rate_NLoS_save_yz minrate(B,pt,q,R_pd,N,p_NLoS_yz_tmp/pt)];
        rate_total_save = [rate_total_save rate_total];
        % total_diff(i,j)=p_LoS/p_NLoS_singlePD;
    end
end
% total_diff = total_diff.';
% thre = 4;
% total_diff_reshape = reshape(total_diff,1,[]);
% LoS_serve_x = [];   LoS_serve_y = [];
% NLoS_serve_x = [];  NLoS_serve_y = [];
% Both_serve_x = [];  Both_serve_y = [];
% for i = 1:length(total_diff_reshape)
%     if total_diff_reshape(i)>thre
%         LoS_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
%         LoS_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
%     elseif total_diff_reshape(i)<1/thre
%         NLoS_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
%         NLoS_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
%     else
%         Both_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
%         Both_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(unique(xd_save), unique(yd_save));
Z_NLoS_power = reshape(p_NLoS_save,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_NLoS_power_xz = reshape(p_NLoS_save_xz,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_NLoS_power_yz = reshape(p_NLoS_save_yz,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_LoS_power = reshape(p_LoS_save,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_power = reshape(p_LoS_save+p_NLoS_save,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_LoS_rate = reshape(rate_LoS_save,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_NLoS_rate_xz = reshape(rate_NLoS_save_xz,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_NLoS_rate_yz = reshape(rate_NLoS_save_yz,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_rate = reshape(rate_total_save,length(yd_save),length(xd_save));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
eval(strcat('X_use_xz = RIS_x_use_'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
eval(strcat('Z_use_xz = RIS_z_use_'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
plot(X_use_xz,Z_use_xz,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
hold on
eval(strcat('X_unuse_xz = RIS_x_unuse_'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
eval(strcat('Z_unuse_xz = RIS_z_unuse_'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
plot(X_unuse_xz,Z_unuse_xz,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
xlim([0,8])
ylim([1,3])
xlabel('X');
ylabel('Z');
title("XZ wall RIS分布图")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
eval(strcat('Y_use_yz = RIS_y_use_yz_'+string(xd_for_plot_yz)+'_'+string(yd_for_plot_yz),';'));
eval(strcat('Z_use_yz = RIS_z_use_yz_'+string(xd_for_plot_yz)+'_'+string(yd_for_plot_yz),';'));
plot(Y_use_yz,Z_use_yz,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
hold on
eval(strcat('Y_unuse_yz = RIS_y_unuse_yz_'+string(xd_for_plot_yz)+'_'+string(yd_for_plot_yz),';'));
eval(strcat('Z_unuse_yz = RIS_z_unuse_yz_'+string(xd_for_plot_yz)+'_'+string(yd_for_plot_yz),';'));
plot(Y_unuse_yz,Z_unuse_yz,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
xlim([0,8])
ylim([1,3])
xlabel('Y');
ylabel('Z');
title("YZ wall RIS分布图")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_power(:)), min(Z_LoS_power(:)));
globalMax = max(max(Z_total_power(:)), max(Z_LoS_power(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_total_power);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('irs分布在相邻两面墙,功率总和');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_LoS_power,3);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_rate(:)), min(Z_LoS_rate(:)));
globalMax = max(max(Z_total_rate(:)), max(Z_LoS_rate(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_total_rate);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('irs分布在相邻两面墙,存在NLoS径时的最低速率');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_LoS_rate,3);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('只有LoS径时的最低速率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_NLoS_power_xz(:)), min(Z_NLoS_power_yz(:)));
globalMax = max(max(Z_NLoS_power_xz(:)), max(Z_NLoS_power_yz(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_NLoS_power_xz);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('xz irs功率');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_NLoS_power_yz);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('yz irs功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_NLoS_rate_xz(:)), min(Z_NLoS_rate_yz(:)));
globalMax = max(max(Z_NLoS_rate_xz(:)), max(Z_NLoS_rate_yz(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_NLoS_rate_xz);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('xz irs的最低速率');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_NLoS_rate_yz);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('yz irs的最低速率');









