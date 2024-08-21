function area = cal_pd_projection_area(point,Mirror_point,angle_ris)
    ris_normal = [sin(angle_ris(1))*cos(angle_ris(2)) sin(angle_ris(1))*sin(angle_ris(2)) cos(angle_ris(1))];
    n1 = ris_normal(1);
    n2 = ris_normal(2);
    n3 = ris_normal(3);
    D_t = -(n1*Mirror_point(1)+n2*Mirror_point(2)+n3*Mirror_point(3));
    point_projection = [];
    for i = 1:size(point,1)
        point_pd = point(i,:);        
        t = (n1*point_pd(1)+n2*point_pd(2)+n3*point_pd(3)+D_t)/(n1^2+n2^2+n3^2);
        point_projection = [point_projection;point_pd(1)-n1*t,point_pd(2)-n2*t,point_pd(3)-n3*t];
    end
    A = point_projection(1,:);
    B = point_projection(2,:);
    C = point_projection(3,:);
    D = point_projection(4,:);
    AD = D-A;
    BC = C-B;
    cos_theta = dot(AD,BC)/(norm(AD)*norm(BC));
    sin_theta = abs(sqrt(1-cos_theta^2));
    area = 0.5*norm(AD)*norm(BC)*sin_theta;
end