function h_gain_NLoS = channel_gain_NLoS_mi(source_point,D_point,reflect_point,source_angle,D_angle,phi_semi,A_PD,f,xi_fov,T,rho_mi)
%%%%%%%% f: internal refractive index       T: gain of the optical filter    
    xs = source_point(1);  ys = source_point(2);  zs = source_point(3);
    xd = D_point(1);  yd = D_point(2);  zd = D_point(3);
    xr = reflect_point(1);  yr = reflect_point(2);  zr = reflect_point(3);

    m = -(log2(cos(phi_semi)))^(-1);

    RD_vectx = xd - xr;     RD_vecty = yd - yr;     RD_vectz = zd - zr;
    RS_vectx = xs - xr;     RS_vecty = ys - yr;     RS_vectz = zs - zr;
    RD_vect = [RD_vectx RD_vecty RD_vectz];
    RS_vect = [RS_vectx RS_vecty RS_vectz];

    RD_norm2 = sqrt(RD_vectx^2+RD_vecty^2+RD_vectz^2);
    RS_norm2 = sqrt(RS_vectx^2+RS_vecty^2+RS_vectz^2);

    N_norm = sqrt(RS_norm2^2+RD_norm2^2+2*(RS_vectx*RD_vectx + RS_vecty*RD_vecty + RS_vectz*RD_vectz));
    N_x = (RS_vectx+RD_vectx)/N_norm;
    N_y = (RS_vecty+RD_vecty)/N_norm;
    N_z = (RS_vectz+RD_vectz)/N_norm;
    % N_x = (RS_vectx+RD_vectx)/sqrt(2+2*(RS_vectx*RD_vectx + RS_vecty*RD_vecty + RS_vectz*RD_vectz));
    % N_y = (RS_vecty+RD_vecty)/sqrt(2+2*(RS_vectx*RD_vectx + RS_vecty*RD_vecty + RS_vectz*RD_vectz));
    % N_z = (RS_vectz+RD_vectz)/sqrt(2+2*(RS_vectx*RD_vectx + RS_vecty*RD_vecty + RS_vectz*RD_vectz));

    source_normal = [sin(source_angle(1))*cos(source_angle(2)) sin(source_angle(1))*sin(source_angle(2)) cos(source_angle(1))];
    D_normal = [sin(D_angle(1))*cos(D_angle(2)) sin(D_angle(1))*sin(D_angle(2)) cos(D_angle(1))];

    phi_RS = acos(dot(-RS_vect,source_normal)/(norm(-RS_vect)*norm(source_normal)));
    theta_RD = acos(dot(RD_vect,-D_normal)/(norm(RD_vect)*norm(-D_normal)));

    G = f^2/(sin(xi_fov)^2);
    A_eff = A_PD *G *T *cos(theta_RD);
    if theta_RD > xi_fov
        h_gain_NLoS = 0;
    else
        h_gain_NLoS = rho_mi*(m+1)*cos(phi_RS)^(m+1)*cos(theta_RD)*A_eff/(2*pi*(RD_norm2+RS_norm2)^2);
    end    
end