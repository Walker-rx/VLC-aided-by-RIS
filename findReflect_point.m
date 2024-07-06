function reflect_point = findReflect_point(AP_point,U_point,Mirror_point,Mirror_angle,wm,hm)
    reflect_point = zeros(3,1);
    % AP_point: AP坐标 [x1, y1, z1]
    % U_point: 接收端坐标 [x2, y2, z2]
    % Mirror_point: RIS中心点 
    % Mirror_angle:[roll yaw]
    % plane_normal: RIS的法向量 [n1, n2, n3]  

    % Counting plane_normal
    roll = Mirror_angle(1);
    yaw = Mirror_angle(2);
    n1 = sin(roll)*cos(yaw);
    n2 = sin(roll)*sin(yaw);
    n3 = cos(roll);

    % 平面方程 Ax+By+Cz+D=0
    D = -(n1*Mirror_point(1)+n2*Mirror_point(2)+n3*Mirror_point(3));

    % 计算对称点AP_prime
    AP_prime = zeros(1,3);
    AP_prime(1) = AP_point(1) - 2*n1*((n1*AP_point(1)+n2*AP_point(2)+n3*AP_point(3)+D)/(n1^2+n2^2+n3^2));
    AP_prime(2) = AP_point(2) - 2*n2*((n1*AP_point(1)+n2*AP_point(2)+n3*AP_point(3)+D)/(n1^2+n2^2+n3^2));
    AP_prime(3) = AP_point(3) - 2*n3*((n1*AP_point(1)+n2*AP_point(2)+n3*AP_point(3)+D)/(n1^2+n2^2+n3^2));

    % 计算与RIS的交点 
    % 线段方程 x=x1+t*(x2-x1)  y=y1+t*(y2-y1)  z=z1+t*(z2-z1)  0<=t<=1
    % AP_prime(x1,y1,z1)  U_point(x2,y2,z2)
    m = U_point(1)-AP_prime(1);
    n = U_point(2)-AP_prime(2);
    p = U_point(3)-AP_prime(3);
    t = -(n1*AP_prime(1)+n2*AP_prime(2)+n3*AP_prime(3)+D)/(n1*m+n2*n+n3*p);
    if (0<t)&&(t<1)
        reflect_point(1) = AP_prime(1) + t*m;
        reflect_point(2) = AP_prime(2) + t*n; 
        reflect_point(3) = AP_prime(3) + t*p;
    else
        reflect_point = [];
    end
    
    % x_limit = abs(wm*cos(yaw-pi/2)/2);
    % y_limit = abs(wm*sin(yaw-pi/2)/2);
    % z_limit = abs(hm*cos(roll-pi/2)/2);
    % 
    % if ~isempty(reflect_point)
    %     if (abs(reflect_point(1)-Mirror_point(1))-x_limit>eps) || (abs(reflect_point(2)-Mirror_point(2))-y_limit>eps) ||(abs(reflect_point(3)-Mirror_point(3))-z_limit>eps)
    %         reflect_point = [];
    %     end
    % end

    if ~isempty(reflect_point)
        if (norm(reflect_point-Mirror_point)-wm/2>eps)
            reflect_point = [];
        end
    end
   
end

