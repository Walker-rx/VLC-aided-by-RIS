function if_intersect = judge_intersect(coor_source,coor_destination,coor_obstacle,angle_obstacle,length_obstacle)  
    % Counting plane_normal
    polar = angle_obstacle(1);
    azi = angle_obstacle(2);
    n1 = sin(polar)*cos(azi);
    n2 = sin(polar)*sin(azi);
    n3 = cos(polar);

    % 平面方程 Ax+By+Cz+D=0
    D = -(n1*coor_obstacle(1)+n2*coor_obstacle(2)+n3*coor_obstacle(3));

    % 计算与obstacle的交点 
    % 线段方程 x=x1+t*(x2-x1)  y=y1+t*(y2-y1)  z=z1+t*(z2-z1)  0<=t<=1
    % coor_source(x1,y1,z1)  coor_destination(x2,y2,z2)
    m = coor_destination(1)-coor_source(1);
    n = coor_destination(2)-coor_source(2);
    p = coor_destination(3)-coor_source(3);
    t = -(n1*coor_source(1)+n2*coor_source(2)+n3*coor_source(3)+D)/(n1*m+n2*n+n3*p);
    if (0<t)&&(t<1)
        intersect_point(1) = coor_source(1) + t*m;
        intersect_point(2) = coor_source(2) + t*n; 
        intersect_point(3) = coor_source(3) + t*p;
    else
        intersect_point = [];
        if_intersect = 0;
        return;
    end
    % z_limit = length_obstacle/2*cos(polar);
    % x_limit = length_obstacle/2*sin(polar)*cos(azi);
    % y_limit = length_obstacle/2*sin(polar)*sin(azi);
    % x_limit = abs(length_obstacle*cos(azi-pi/2)/2);
    % y_limit = abs(length_obstacle*sin(azi-pi/2)/2);
    % z_limit = abs(length_obstacle*cos(polar-pi/2)/2);
    % 
    % if (abs(intersect_point(1)-coor_obstacle(1))-x_limit)>eps && (abs(intersect_point(2)-coor_obstacle(2))-y_limit)>eps && (abs(intersect_point(3)-coor_obstacle(3))-z_limit)>eps
    %     if_intersect = 0;
    % else
    %     if_intersect = 1;
    % end
    if (norm(intersect_point-coor_obstacle)-length_obstacle/2>eps)
        if_intersect = 0;
    else
        if_intersect = 1;
    end
end