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

hd = 3;
pd_polar = pi; pd_zai = 0;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
angle_PD = [pd_polar,pd_zai];

wm = 4/100; hm = 4/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
xr = 1;   zr = 1.5; yr = 1;
total_length = 7-xr;  total_hight = 2.5-zr;
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
xd_for_plot_xz_oppostie = 9;  yd_for_plot_xz_oppostie = 9;
p_LoS_save = [];
p_NLoS_save = [];
rate_LoS_save = [];
rate_total_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

% coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(xd_save)
    xd = xd_save(i);
    for j = 1:length(yd_save)
        str_diff=strcat('diff_'+string(i)+'_'+string(j));
        eval(strcat(str_diff+'=[];'));
        str_x_use=strcat('RIS_x_use_'+string(i)+'_'+string(j));
        eval(strcat(str_x_use+'=[];'));
        str_z_use=strcat('RIS_z_use_'+string(i)+'_'+string(j));
        eval(strcat(str_z_use+'=[];'));
        str_x_unuse=strcat('RIS_x_unuse_'+string(i)+'_'+string(j));
        eval(strcat(str_x_unuse+'=[];'));
        str_z_unuse=strcat('RIS_z_unuse_'+string(i)+'_'+string(j));
        eval(strcat(str_z_unuse+'=[];'));
        
        str_diff_xz_oppostie=strcat('diff_xz_oppostie_'+string(i)+'_'+string(j));
        eval(strcat(str_diff_xz_oppostie+'=[];'));
        str_x_use_xz_oppostie=strcat('RIS_x_use_xz_oppostie_'+string(i)+'_'+string(j));
        eval(strcat(str_x_use_xz_oppostie+'=[];'));
        str_z_use_xz_oppostie=strcat('RIS_z_use_xz_oppostie_'+string(i)+'_'+string(j));
        eval(strcat(str_z_use_xz_oppostie+'=[];'));
        str_x_unuse_xz_oppostie=strcat('RIS_x_unuse_xz_oppostie_'+string(i)+'_'+string(j));
        eval(strcat(str_x_unuse_xz_oppostie+'=[];'));
        str_z_unuse_xz_oppostie=strcat('RIS_z_unuse_xz_oppostie_'+string(i)+'_'+string(j));
        eval(strcat(str_z_unuse_xz_oppostie+'=[];'));
        
        yd = yd_save(j);
        coor_PD = [xd,yd,hd].';
        p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);        
        p_NLoS_singlePD = 0;
        rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);

                coor_RIS_xz_oppostie = [xr+wm/2+(jj-1)*(wm+delta_wm), 8, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS_xz_oppostie = findReflect_angle_xz_oppostie(coor_led,coor_PD,coor_RIS_xz_oppostie);
                if all((angle_RIS > 0) & (angle_RIS < pi))
                    coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                    if all(abs(coor_reflect-coor_RIS)<1e-7)
                        p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        diff_tmp = p_LoS/p_NLoS;
                        eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;
                        if diff_tmp == Inf
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%% _xz_oppostie %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (angle_RIS_xz_oppostie(1) > 0) && (angle_RIS_xz_oppostie(1) < pi) && (angle_RIS_xz_oppostie(2) > -pi) && (angle_RIS_xz_oppostie(2) < 0)
                    coor_reflect_xz_oppostie = findReflect_point(coor_led,coor_PD,coor_RIS_xz_oppostie,angle_RIS_xz_oppostie,wm,hm);
                    if all(abs(coor_reflect_xz_oppostie-coor_RIS_xz_oppostie)<1e-7)
                        p_NLoS_xz_oppostie = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect_xz_oppostie,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        diff_tmp_xz_oppostie = p_LoS/p_NLoS_xz_oppostie;
                        eval(strcat(str_diff_xz_oppostie+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp_xz_oppostie)+';'));
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS_xz_oppostie;
                        if diff_tmp_xz_oppostie == Inf
                            eval(strcat(str_x_unuse_xz_oppostie,'=[',str_x_unuse_xz_oppostie,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_unuse_xz_oppostie,'=[',str_z_unuse_xz_oppostie,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                        else
                            eval(strcat(str_x_use_xz_oppostie,'=[',str_x_use_xz_oppostie,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                            eval(strcat(str_z_use_xz_oppostie,'=[',str_z_use_xz_oppostie,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));                            
                        end                       
                    else
                        diff_tmp_xz_oppostie = NaN;
                        fprintf("NaN happen in xz oppostie wall, i=%d,j=%d,ii=%d,jj=%d \n",i,j,ii,jj);
                        eval(strcat(str_diff_xz_oppostie+'('+string(ii)+','+string(jj)+')'+'=NaN;'));
                        eval(strcat(str_x_unuse_xz_oppostie,'=[',str_x_unuse_xz_oppostie,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                        eval(strcat(str_z_unuse_xz_oppostie,'=[',str_z_unuse_xz_oppostie,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                    end                    
                else
                    diff_tmp_xz_oppostie = NaN;
                    fprintf("NaN happen in yz wall, i=%d,j=%d,ii=%d,jj=%d \n",i,j,ii,jj);
                    eval(strcat(str_diff_xz_oppostie+'('+string(ii)+','+string(jj)+')'+'=NaN;'));
                    eval(strcat(str_x_unuse_xz_oppostie,'=[',str_x_unuse_xz_oppostie,','+string(xr+wm/2+(jj-1)*(wm+delta_wm))+'];'));
                    eval(strcat(str_z_unuse_xz_oppostie,'=[',str_z_unuse_xz_oppostie,','+string(zr+hm/2+(ii-1)*(hm+delta_hm))+'];'));
                end

                
            end
        end

        rate_total = minrate(B,pt,q,R_pd,N,(p_LoS+p_NLoS_singlePD)/pt);
        p_LoS_save = [p_LoS_save p_LoS];
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
        rate_LoS_save = [rate_LoS_save rate_LoS];
        rate_total_save = [rate_total_save rate_total];
        total_diff(i,j)=p_LoS/p_NLoS_singlePD;
    end
end
total_diff = total_diff.';
thre = 4;
total_diff_reshape = reshape(total_diff,1,[]);
LoS_serve_x = [];   LoS_serve_y = [];
NLoS_serve_x = [];  NLoS_serve_y = [];
Both_serve_x = [];  Both_serve_y = [];
for i = 1:length(total_diff_reshape)
    if total_diff_reshape(i)>thre
        LoS_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
        LoS_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
    elseif total_diff_reshape(i)<1/thre
        NLoS_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
        NLoS_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
    else
        Both_serve_x(end+1) = xd_save(floor((i-1)/length(yd_save))+1);
        Both_serve_y(end+1) = yd_save(mod(i-1,length(xd_save))+1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X, Y] = meshgrid(unique(xd_save), unique(yd_save));
Z_NLoS_power = reshape(p_NLoS_save,length(yd_save),length(xd_save));
% figure;
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
% contourf(X, Y, Z_NLoS_power, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('NLOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_LoS_power = reshape(p_LoS_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS_power(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('LOS功率');
% title('LOS功率');

% figure
% contourf(X, Y, Z_LoS_power, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_power = reshape(p_LoS_save+p_NLoS_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_power(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('功率总和');
% title('功率总和');

% figure
% contourf(X, Y, Z_total_power, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('功率总和');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_diff = total_diff;
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_diff(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('LOS功率/NLOS功率');
% title('LOS功率/NLOS功率');

figure
contourf(X, Y, Z_diff, 6);
colorbar;

xlabel('X');
ylabel('Y');
title('LOS功率/NLOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_diff_inverse = 1./total_diff;
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_diff_inverse(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('NLOS功率/LOS功率');
% title('NLOS功率/LOS功率');

% figure
% contourf(X, Y, Z_diff_inverse, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('NLOS功率/LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_diff = reshape((p_LoS_save+p_NLoS_save)./p_LoS_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_diff(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('功率总和/LOS功率');
% title('功率总和/LOS功率');

% figure
% contourf(X, Y, Z_total_diff, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('功率总和/LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_LoS_rate = reshape(rate_LoS_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS_rate(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('最低速率');
% title('只有LoS径时的最低速率');

% figure
% contourf(X, Y, Z_LoS_rate, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('只有LoS径时的最低速率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z_total_rate = reshape(rate_total_save,length(yd_save),length(xd_save));
% figure;
% trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_rate(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
% colorbar;  % 添加颜色条
% 
% % 设置坐标轴标签
% xlabel('X');
% ylabel('Y');
% zlabel('最低速率');
% title('存在NLoS径时的最低速率');

% figure
% contourf(X, Y, Z_total_rate, 6);
% colorbar;
% 
% xlabel('X');
% ylabel('Y');
% title('存在NLoS径时的最低速率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
if ~isempty(LoS_serve_x)
    plot(LoS_serve_x,LoS_serve_y,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
    hold on
    plot(Both_serve_x,Both_serve_y,'go', 'MarkerFaceColor', 'green', 'MarkerSize', 8)
    if ~isempty(NLoS_serve_x)
        hold on
        plot(NLoS_serve_x,NLoS_serve_y,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
        legend('LoS Service', 'Both Services', 'NLoS Service');
    else
        legend('LoS Service', 'Both Services');
    end
else
    if ~isempty(NLoS_serve_x)
        plot(NLoS_serve_x,NLoS_serve_y,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
        hold on
        plot(Both_serve_x,Both_serve_y,'go', 'MarkerFaceColor', 'green', 'MarkerSize', 8)
        legend('NLoS Service', 'Both Services');
    else
        plot(Both_serve_x,Both_serve_y,'go', 'MarkerFaceColor', 'green', 'MarkerSize', 8)
        legend('Both Services');
    end
end


xlim([0,8])
ylim([0,8])
xlabel('X');
ylabel('Y');
title("服务分布图")
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
eval(strcat('Y_use_xz_oppostie = RIS_x_use_xz_oppostie_'+string(xd_for_plot_xz_oppostie)+'_'+string(yd_for_plot_xz_oppostie),';'));
eval(strcat('Z_use_xz_oppostie = RIS_z_use_xz_oppostie_'+string(xd_for_plot_xz_oppostie)+'_'+string(yd_for_plot_xz_oppostie),';'));
plot(Y_use_xz_oppostie,Z_use_xz_oppostie,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
hold on
eval(strcat('Y_unuse_xz_oppostie = RIS_x_unuse_xz_oppostie_'+string(xd_for_plot_xz_oppostie)+'_'+string(yd_for_plot_xz_oppostie),';'));
eval(strcat('Z_unuse_xz_oppostie = RIS_z_unuse_xz_oppostie_'+string(xd_for_plot_xz_oppostie)+'_'+string(yd_for_plot_xz_oppostie),';'));
plot(Y_unuse_xz_oppostie,Z_unuse_xz_oppostie,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
xlim([0,8])
ylim([1,3])
xlabel('Y');
ylabel('Z');
title("XZ oppostie wall RIS分布图")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_power(:)), min(Z_LoS_power(:)));
globalMax = max(max(Z_total_power(:)), max(Z_LoS_power(:)));

subplot(2, 1, 1);
contourf(X, Y, Z_total_power);
caxis([globalMin globalMax]);
colorbar('location', 'EastOutside');
xlabel('X');
ylabel('Y');
title('功率总和');

subplot(2, 1, 2);
contourf(X, Y, Z_LoS_power);
caxis([globalMin globalMax]);
colorbar('location', 'EastOutside');
xlabel('X');
ylabel('Y');
title('LOS功率');
colormap(jet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_power(:)), min(Z_LoS_power(:)));
globalMax = max(max(Z_total_power(:)), max(Z_LoS_power(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_total_power);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('irs分布在对面两面墙,功率总和');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_LoS_power);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_rate(:)), min(Z_LoS_rate(:)));
globalMax = max(max(Z_total_rate(:)), max(Z_LoS_rate(:)));

subplot(2, 1, 1);
contourf(X, Y, Z_total_rate);
caxis([globalMin globalMax]);
colorbar('location', 'EastOutside');
xlabel('X');
ylabel('Y');
title('存在NLoS径时的最低速率');

subplot(2, 1, 2);
contourf(X, Y, Z_LoS_rate);
caxis([globalMin globalMax]);
colorbar('location', 'EastOutside');
xlabel('X');
ylabel('Y');
title('只有LoS径时的最低速率');
colormap(jet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
globalMin = min(min(Z_total_rate(:)), min(Z_LoS_rate(:)));
globalMax = max(max(Z_total_rate(:)), max(Z_LoS_rate(:)));

subplot(2, 1, 1);
[C,h]=contour(X, Y, Z_total_rate);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('irs分布在对面两面墙,存在NLoS径时的最低速率');

subplot(2, 1, 2);
[C,h]=contour(X, Y, Z_LoS_rate);
clabel(C,h)
xlabel('X');
ylabel('Y');
title('只有LoS径时的最低速率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS_rate(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
hold on
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_rate(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
xlabel('X');
ylabel('Y');
title('IRS对最低速率的提升');
colormap(jet);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_LoS_power(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
hold on
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z_total_power(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
xlabel('X');
ylabel('Y');
title('IRS对功率的提升');
colormap(jet);








