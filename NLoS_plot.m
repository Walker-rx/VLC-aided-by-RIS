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

% k = 30; l = 1;

xd_save = 1:0.5:7;
yd_save = 1:0.5:7;
p_LoS_save = [];
p_NLoS_save = [];
total_diff = zeros(length(xd_save),length(yd_save));

xr = 1.5;zr = 2;
k = 12; l = 90;

coor_RIS = [xr+wm/2+(l-1)*(wm+delta_wm), 0, zr+hm/2+(k-1)*(hm+delta_hm)].';
for i = 1:length(xd_save)
    xd = xd_save(i);
    for j = 1:length(yd_save)
        str1=strcat('diff_'+string(i)+'_'+string(j));
        eval(strcat(str1+'=[];'));
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
                        eval(strcat(str1+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                        p_NLoS_singlePD = p_NLoS_singlePD + p_NLoS;
                    else
                        diff_tmp = NaN;
                        fprintf("NaN happen in i=%d,j=%d,ii=%d,jj=%d",i,j,ii,jj);
                        eval(strcat(str1+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)+';'));
                    end                    
                else
                    diff_tmp = NaN;
                    eval(strcat(str1+'('+string(ii)+','+string(jj)+')'+'='+string(diff_tmp)));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = reshape(p_LoS_save,length(yd_save),length(xd_save));
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('LOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = total_diff;
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('LOS功率/NLOS功率');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = 1./total_diff;
figure;
trisurf(delaunay(X(:), Y(:)), X(:), Y(:), Z(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
colorbar;  % 添加颜色条

% 设置坐标轴标签
xlabel('X');
ylabel('Y');
zlabel('NLOS功率/LOS功率');








