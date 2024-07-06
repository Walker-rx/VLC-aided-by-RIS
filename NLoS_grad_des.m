clear
close all

syms alpha_RIS beta_RIS % alpha:polar angle ; beta: azimuth angle
syms ys
syms xr yr zr
syms xd yd hd
syms alpha_source beta_source alpha_PD beta_PD
syms f xi_fov A_PD T phi_semi rho_mi
pt = 20;

coor_led = [0,ys,0].';
coor_RIS = [xr yr zr].';
coor_PD = [xd,yd,hd].';

%%%%%%%%%  Find reflect point %%%%%%%%%%%
N_hat = [sin(alpha_RIS)*cos(beta_RIS) sin(alpha_RIS)*sin(beta_RIS) cos(alpha_RIS)].';
vect_RS = coor_led-coor_RIS;
coor_Sprime = coor_led - 2*N_hat.'*vect_RS*N_hat;
vect_RSprime = coor_Sprime - coor_RIS;

vect_SprimeP = coor_PD-coor_Sprime;
t = -1*(N_hat.'*vect_RSprime)/(N_hat.'*vect_SprimeP);
coor_reflect = coor_Sprime+t*vect_SprimeP;

%%%%%%%%% Cal NLoS gain %%%%%%%%%%%%%%%%%
vect_reflect_P = coor_PD - coor_reflect;
vect_reflect_S = coor_led - coor_reflect;
normal_source = [sin(alpha_source)*cos(beta_source) sin(alpha_source)*sin(beta_source) cos(alpha_source)];
normal_PD = [sin(alpha_PD)*cos(beta_PD) sin(alpha_PD)*sin(beta_PD) cos(alpha_PD)];

norm2_vect_reflectP = norm(vect_reflect_P);
norm2_vect_reflectS = norm(vect_reflect_S);

phi_RS = acos(dot(-vect_reflect_S,normal_source)/(norm(-vect_reflect_S)*norm(normal_source)));
theta_RD = acos(dot(vect_reflect_P,-normal_PD)/(norm(vect_reflect_P)*norm(-normal_PD)));

G = f^2/(sin(xi_fov)^2);
A_eff = A_PD *G *T *cos(theta_RD);
m = -(log2(cos(phi_semi)))^(-1);

h_gain_NLoS = rho_mi*(m+1)*cos(phi_RS)^(m+1)*cos(theta_RD)*A_eff/(2*pi*(norm2_vect_reflectP+norm2_vect_reflectS)^2);
p_NLoS = pt*h_gain_NLoS;
p_NLoS = matlabFunction(p_NLoS);
diff_alpha = diff(p_NLoS,alpha_RIS);
diff_alpha = matlabFunction(diff_alpha);
diff_beta = diff(p_NLoS,beta_RIS);
diff_beta = matlabFunction(diff_beta);

[polar_point,azi_point,P_max] = gradient_descent(p_NLoS,diff_alpha,diff_beta ,...
    1.000000000000000e-06, 1, pi, pi*0.4953, 0, 0, 0.4970*pi, 0, 1.500000000000000, 3, pi/3 ,0.8 ,0, 85/180*pi, 0, 2, 0, 2, 1.5);
