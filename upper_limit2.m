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

xd = 3; yd = 3; hd = 3;
pd_polar = pi; pd_zai = 0;
wd = 1/100; ld = 1/100;
A_PD = wd*ld;
angle_PD = [pd_polar,pd_zai];
coor_PD = [xd,yd,hd].';


% xr = 1;   zr = 1.5;
% total_length = 7-xr;  total_hight = 2.5-zr;
% xr = 2.5;   zr = 1.5;
% total_length = 3.5-xr;  total_hight = 2.5-zr;

delta_wm = 0; delta_hm = 0;
total_length_loop = 1:0.5:6;
zr = 1.5; total_hight = 2.5-zr;

area_save = [];
total_area_save = [];
num_save = [];
p_NLoS_save = [];
rate_NLoS_save = [];

p_LoS = pt*channel_gain_LoS(coor_led,coor_PD,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T); 
rate_LoS = minrate(B,pt,q,R_pd,N,p_LoS/pt);
k=300; l=200;
for iii = 1:length(total_length_loop)
    total_length = total_length_loop(iii);
    xr = 3-0.5*total_length;
    hm = total_hight/k;
    wm = total_length/l;
    area = hm*wm;       
    num = k*l;
    p_NLoS_singlePD = 0;  
    for i = 1:k
        for j = 1:l
            coor_RIS = [xr+wm/2+(j-1)*(wm+delta_wm), 0, zr+hm/2+(i-1)*(hm+delta_hm)].';
            angle_RIS = findReflect_angle_xz(coor_led,coor_PD,coor_RIS);
            if all((angle_RIS > 0) & (angle_RIS < pi))
                coor_reflect = findReflect_point(coor_led,coor_PD,coor_RIS,angle_RIS,wm,hm);
                pd_point = [xd-wd/2,yd-ld/2,hd;  xd-wd/2,yd+ld/2,hd;  xd+wd/2,yd-ld/2,hd; xd+wd/2,yd+ld/2,hd;];
                % pd_area = cal_pd_projection_area(pd_point,coor_RIS,angle_RIS);
                if all(abs(coor_reflect-coor_RIS)<1e-7)
                    if area>1e-4
                        p_NLoS = pt*channel_gain_NLoS_mi(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,phi_semi,A_PD,f,xi_fov,T,rho_mi);
                    else
                        p_NLoS = pt*channel_gain_NLoS_mi_upper(coor_led,coor_PD,coor_reflect,angle_led,angle_PD,angle_RIS,phi_semi,area,f,xi_fov,T,rho_mi);
                    end                     
                    p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;                                         
                else                      
                    fprintf("reflect error happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);                        
                end                    
            else
                fprintf("angle error happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
            end
        end
    end
    rate_NLoS = minrate(B,pt,q,R_pd,N,p_NLoS_singlePD/pt);  
    p_NLoS_save = [p_NLoS_save,p_NLoS_singlePD];
    rate_NLoS_save = [rate_NLoS_save,rate_NLoS];
end

figure
plot(total_length_loop,p_NLoS_save)
title('NLoS径功率上限随面积变化')
figure
plot(total_length_loop,rate_NLoS_save)
title('NLoS径速率上限随面积变化')


