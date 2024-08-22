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
xr_max_loop = 1:0.5:8;   zr_max_loop = 0.5:0.2:3;
total_length = 1;  total_hight = 0.5;
max_k = total_hight/hm;    max_l=total_length/wm;

k = 8; l = 20;
% k = 16; l = 40;
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

xd_save = 1:0.5:7;
yd_save = 1:0.5:7;
xd_for_plot = 1;  zd_for_plot = 5;
p_LoS_save = [];
p_NLoS_save = [];
rate_LoS_save = [];
rate_total_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

coor_PD = [7;3;3];
p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T);
rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);
% coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(xr_max_loop)
    xr = xr_max_loop(i)-total_length;
    for j = 1:length(zr_max_loop)
        zr = zr_max_loop(j)-total_hight;    
                
        p_NLoS_singlePD = 0;
        
        for ii = 1:k
            for jj = 1:l
                coor_RIS = [xr+wm/2+(jj-1)*(wm+delta_wm), 0, zr+hm/2+(ii-1)*(hm+delta_hm)].';
                angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);
                if all((angle_RIS > 0) & (angle_RIS < pi))
                    coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                    if all(abs(coor_reflect-coor_RIS)<1e-7)
                        p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);                   
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;
                                          
                    else
                        diff_tmp = NaN;
                        fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);

                    end                    
                else
                    diff_tmp = NaN;
                    fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d \n",i,j,ii,jj);
                   
                end
            end
        end

        rate_total = minrate(B,pt,q,R_pd,N,(p_LoS+p_NLoS_singlePD)/pt);
        p_NLoS_save = [p_NLoS_save p_NLoS_singlePD];
        rate_total_save = [rate_total_save rate_total];
    end
end
[X, Y] = meshgrid(unique(xr_max_loop-total_length/2), unique(zr_max_loop-total_hight/2));
Z_total_rate = reshape(rate_total_save,length(zr_max_loop),length(xr_max_loop));
[C,h]=contour(X, Y, Z_total_rate);
clabel(C,h);
xlabel('IRS部署区域中心点的X坐标');
ylabel('IRS部署区域中心点的Z坐标');
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
title('PD从不同区域IRS处接收的速率下界(bit/s)')