function reflect_angle = findReflect_angle_yz(LED_point,PD_point,Mirror_point)
    RS_vect = LED_point - Mirror_point;
    RD_vect = PD_point - Mirror_point;
    RS_normal = RS_vect./norm(RS_vect);
    RD_normal = RD_vect./norm(RD_vect);
    % n = RS_vect + RD_vect;
    aa = RS_normal+RD_normal;
    n = [1,aa(2)/aa(1),aa(3)/aa(1)].';
    
    n_norm = norm(n);
    n_normal = n./n_norm;
    
    sigma1 = acos(dot(RS_vect,RD_vect)/(norm(RS_vect)*norm(RD_vect)));
    sigma2 = acos(dot(n_normal,RD_vect)/(norm(n_normal)*norm(RD_vect)));
    sigma3 = acos(dot(RS_vect, n_normal )/(norm(RS_vect)*norm(n_normal )));
    if abs(sigma1-2*sigma2) >1e-13 || abs(sigma2-sigma3) >1e-13 
        n_normal = -1*n_normal;
    end
    sigma2 = acos(dot(n_normal,RD_vect)/(norm(n_normal)*norm(RD_vect)));
    sigma3 = acos(dot(RS_vect, n_normal )/(norm(RS_vect)*norm(n_normal )));
    
    cos_polar =  n_normal(3);
    tan_azi = n_normal(2)/n_normal(1);
    
    angle_polar = acos(cos_polar);
    angle_azi = atan(tan_azi);
    reflect_angle = [angle_polar,angle_azi];
end