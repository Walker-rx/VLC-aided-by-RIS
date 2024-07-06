function h_gain_LoS = channel_gain_LoS(source_point,PD_point,source_angle,PD_angle,phi_semi,A_PD,f,xi_fov,T)
%%%%%%%% source_angle PD_angle:(polar angle , azimuth angle)
%%%%%%%% (xa, ya, za): AP       (xu, yu, zu): user     (xk, yk, zk): Mirror      (xn, yn, zn): LC   
%%%%%%%% phi_semi: semi-angle at half power of the LED        A_PD: detector physical area of the photo-diode
%%%%%%%% xi: the incident angle       phi: the radiation angle      xi_fov: FoV of the PD
%%%%%%%% f: internal refractive index       T: gain of the optical filter  
%%%%%%%% alpha: receiver's polar angle         beta: receiver's azimuth angle
    xa = source_point(1);  ya = source_point(2);  za = source_point(3);
    xu = PD_point(1);  yu = PD_point(2);  zu = PD_point(3);
    % xk = Mirror_point(1);  yk = Mirror_point(2);  zk = Mirror_point(3);
    alpha = PD_angle(1);  beta = PD_angle(2);
    
    LoS_normal = [xu-xa yu-ya zu-za];
    AP_normal = [sin(source_angle(1))*cos(source_angle(2)) sin(source_angle(1))*sin(source_angle(2)) cos(source_angle(1))];
    U_normal = [sin(alpha)*cos(beta) sin(alpha)*sin(beta) cos(alpha)];

    xi = acos(dot(LoS_normal,-U_normal)/(norm(LoS_normal)*norm(-U_normal)));
    phi = acos(dot(LoS_normal,AP_normal)/(norm(LoS_normal)*norm(AP_normal)));
    

%%%%%%%%%%%%%%%%%%%%%%    LoS    %%%%%%%%%%%%%%%%%%%%%%   
    m = -(log2(cos(phi_semi)))^(-1);
    d_LoS = sqrt((xa-xu)^2+(ya-yu)^2+(za-zu)^2);
    G = f^2/(sin(xi_fov)^2);
    cos_xi = ((xa-xu)/d_LoS)*cos(beta)*sin(alpha) + ((ya-yu)/d_LoS)*sin(beta)*sin(alpha) + ((za-zu)/d_LoS)*cos(alpha);
    if xi > xi_fov
        h_gain_LoS = 0;
    else
        h_gain_LoS = (m+1) *A_PD *G *T *cos(phi)^m *cos_xi /(2*pi*d_LoS^2);
    end

% %%%%%%%%%%%%%%%%%%%%%%    Mirror NLoS    %%%%%%%%%%%%%%%%%%%%%%
% %%%%%%% dk_a: link distance for AP-Mirror      dk_u:  link distance for Mirror-user
% %%%%%%% phi_ka: radiation angle for the path from the AP to Mirror       xi_ka: incident angle for the path from AP to Mirror
% %%%%%%% phi_uk: angle of irradiance for the path from Mirror towards user u
% %%%%%%% xi_uk: incident angle of the reflected signal from Mirror to user u
% %%%%%%% gamma:  yaw angles of Mirror            omega: roll angles of Mirror
% %%%%%%% 不止cos(phi),cos(xi)也可以用坐标表示？
%     gamma = Mirror_angle(2);  omega = Mirror_angle(1);
%     incident_normal = [xa-reflect_point(1) ya-reflect_point(2) za-reflect_point(3)];
%     emerging_normal = [xu-reflect_point(1) yu-reflect_point(2) zu-reflect_point(3)];
%     Mirror_normal = [cos(omega)*sin(gamma) cos(omega)*cos(gamma) sin(omega)];   
%     phi_ka = acos(dot(-incident_normal,AP_normal)/(norm(-incident_normal)*norm(AP_normal)));
%     xi_ka = acos(dot(incident_normal,Mirror_normal)/(norm(incident_normal)*norm(Mirror_normal)));
%     phi_uk = acos(dot(emerging_normal,Mirror_normal)/(norm(emerging_normal)*norm(Mirror_normal)));
%     xi_uk = acos(dot(-emerging_normal,U_normal)/(norm(-emerging_normal)*norm(U_normal)));
% 
%     dk_a = sqrt((xa-xk)^2+(ya-yk)^2+(za-zk)^2);
%     du_k = sqrt((xk-xu)^2+(yk-yu)^2+(zk-zu)^2);
%     cos_phi_uk = ((xk-xu)/du_k)*sin(gamma)*cos(omega) +((yk-yu)/du_k)*cos(gamma)*cos(omega) +((zk-zu)/du_k)*sin(omega);
%     if xi_uk> xi_fov
%         H_NLoS_Mirror = 0;
%     else
%         H_NLoS_Mirror = (m+1) *A_PD *Ak *G *T *cos(phi_ka)^m *cos(xi_ka) *cos_phi_uk *cos(xi_uk) /(2*pi^2*dk_a^2*du_k^2);
%     end
%     H_NLoS_Mirror = H_NLoS_Mirror*rho;
% 
%     h_gain_LoS = H_Los;

end


