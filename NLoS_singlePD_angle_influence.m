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
xr = 1;   zr = 1.5;
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
xd_for_count = 8;  yd_for_count = 9;
hd = 3;
coor_PD = [xd_save(xd_for_count),yd_save(yd_for_count),hd].';

% pd_polar_min = 0; pd_polar_max = pi;
% pd_polar_save = pd_polar_min : (pd_polar_max-pd_polar_min)/40 : pd_polar_max; 
pd_polar_save = 3*pi/4; 
pd_azi_min = 0; pd_azi_max = 2*pi;
pd_azi_save = pd_azi_min : (pd_azi_max-pd_azi_min)/40 : pd_azi_max;
polar_for_plot = 10;  azi_for_plot = 10;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;

p_LoS_save = [];
p_NLoS_save = [];
rate_LoS_save = [];
rate_total_save = [];

% coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(pd_polar_save)
    pd_polar = pd_polar_save(i);
    for j = 1:length(pd_azi_save)
        pd_azi = pd_azi_save(j);
        angle_PD = [pd_polar,pd_azi];

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
        
        p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);        
        p_NLoS_singlePD = 0;
        rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);
                if all((angle_RIS > 0) & (angle_RIS < pi))
                    coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                    if all(abs(coor_reflect-coor_RIS)<1e-7)
                        p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                        % diff_tmp = p_LoS/p_NLoS;
                        % eval(strcat(str_diff+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
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
        p_LoS_save = [p_LoS_save p_LoS];
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
        rate_LoS_save = [rate_LoS_save rate_LoS];
        rate_total_save = [rate_total_save rate_total];

    end
end
LoS_plot = figure;
PolarPlot_3d(LoS_plot,pd_polar_save,pd_azi_save,p_LoS_save,coor_PD,'Polar','Azimuth','P_{LoS}','LoS径功率随角度变化')
ax1=gca;
LoS_plot_zlim=ax1.ZLim;
NLoS_plot = figure;
PolarPlot_3d(NLoS_plot,pd_polar_save,pd_azi_save,p_NLoS_save+p_LoS_save,coor_PD,'Polar','Azimuth','P_{NLoS}','NLoS径功率随角度变化')
ax2=gca;
NLoS_plot_zlim=ax2.ZLim;
z_min = min(LoS_plot_zlim(1),NLoS_plot_zlim(1));
z_max = max(LoS_plot_zlim(2),NLoS_plot_zlim(2));
figure(LoS_plot)
zlim([z_min,z_max])
figure(NLoS_plot)
zlim([z_min,z_max])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% eval(strcat('X_use = RIS_x_use_'+string(polar_for_plot)+'_'+string(azi_for_plot),';'));
% eval(strcat('Z_use = RIS_z_use_'+string(polar_for_plot)+'_'+string(azi_for_plot),';'));
% plot(X_use,Z_use,'ro', 'MarkerFaceColor', 'red', 'MarkerSize', 8)
% hold on
% eval(strcat('X_unuse = RIS_x_unuse_'+string(polar_for_plot)+'_'+string(azi_for_plot),';'));
% eval(strcat('Z_unuse = RIS_z_unuse_'+string(polar_for_plot)+'_'+string(azi_for_plot),';'));
% plot(X_unuse,Z_unuse,'bo', 'MarkerFaceColor', 'blue', 'MarkerSize', 8)
% xlim([0,8])
% ylim([1,3])
% xlabel('X');
% ylabel('Z');
% title("XZ wall RIS分布图")










