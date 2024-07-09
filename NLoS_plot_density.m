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

wm = 4/100; hm = 4/100; delta_wm = 1/2*wm; delta_hm = 1/2*hm;
xr = 1;   zr = 1.5;
total_length = 7-xr;  total_hight = 2.5-zr;
max_k = total_hight/hm;    max_l=total_length/wm;

k = 5; l = 20;
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
p_LoS_save = [];
p_NLoS_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

% coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(xd_save)
    xd = xd_save(i);
    for j = 1:length(yd_save)
        str_diff=strcat('diff_'+string(i)+'_'+string(j));
        eval(strcat(str_diff+'=[];'));
        str_x_use=strcat('RIS_x_use'+string(i)+'_'+string(j));
        eval(strcat(str_x_use+'=[];'));
        str_z_use=strcat('RIS_z_use'+string(i)+'_'+string(j));
        eval(strcat(str_z_use+'=[];'));
        str_x_unuse=strcat('RIS_x_unuse'+string(i)+'_'+string(j));
        eval(strcat(str_x_unuse+'=[];'));
        str_z_unuse=strcat('RIS_z_unuse'+string(i)+'_'+string(j));
        eval(strcat(str_z_unuse+'=[];'));
        yd = yd_save(j);
        coor_PD = [xd,yd,hd].';
        p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);        
        p_NLoS_singlePD = 0;
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle(coor_led,coor_PD,coor_RIS);
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
            end
        end
        p_LoS_save = [p_LoS_save p_LoS];
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
        total_diff(i,j)=p_LoS/p_NLoS_singlePD;
    end
end
[X, Y] = meshgrid(unique(xd_save), unique(yd_save));
Z = reshape(p_NLoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('NLOS功率');
title('NLOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = reshape(p_LoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('LOS功率');
title('LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = reshape(p_LoS_save+p_NLoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('功率总和');
title('功率总和');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = total_diff;
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('LOS功率/NLOS功率');
title('LOS功率/NLOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = 1./total_diff;
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('NLOS功率/LOS功率');
title('NLOS功率/LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = reshape((p_LoS_save+p_NLoS_save)./p_LoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('功率总和/LOS功率');
title('功率总和/LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
eval(strcat('X = RIS_x_use'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
eval(strcat('Z = RIS_z_use'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
plot(X,Z,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
hold on
eval(strcat('X = RIS_x_unuse'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
eval(strcat('Z = RIS_z_unuse'+string(xd_for_plot)+'_'+string(yd_for_plot),';'));
plot(X,Z,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
xlim([0,8])
ylim([1,3])
xlabel('X');
ylabel('Z');
title("RIS分布图")








