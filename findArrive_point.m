function arrive_point = findArrive_point(AP_point, Mirror_point, Mirror_angle, z0)  
    arrive_point = zeros(3,1);
    roll = Mirror_angle(1);
    yaw = Mirror_angle(2);
    % 计算镜面法向量
    dx = sin(roll) * cos(yaw);
    dy = sin(roll) * sin(yaw);
    dz = cos(roll);
    
    % 计算镜面单位法向量
    unit_dx = dx / sqrt(dx^2 + dy^2 + dz^2);
    unit_dy = dy / sqrt(dx^2 + dy^2 + dz^2);
    unit_dz = dz / sqrt(dx^2 + dy^2 + dz^2);
    unit_RIS = [unit_dx,unit_dy,unit_dz].';
    
    % 计算入射光线方向向量
    vector_to_mirror_center = [Mirror_point(1) - AP_point(1); Mirror_point(2) - AP_point(2); Mirror_point(3) - AP_point(3)];
    
    % 计算反射光线的方向向量
    reflect_dx = vector_to_mirror_center(1) - 2 * unit_dx * dot(vector_to_mirror_center, unit_RIS);
    reflect_dy = vector_to_mirror_center(2) - 2 * unit_dy * dot(vector_to_mirror_center, unit_RIS);
    reflect_dz = vector_to_mirror_center(3) - 2 * unit_dz * dot(vector_to_mirror_center, unit_RIS);
    
    % 计算反射光线的单位方向向量
    unit_reflect_dx = reflect_dx / sqrt(reflect_dx^2 + reflect_dy^2 + reflect_dz^2);
    unit_reflect_dy = reflect_dy / sqrt(reflect_dx^2 + reflect_dy^2 + reflect_dz^2);
    unit_reflect_dz = reflect_dz / sqrt(reflect_dx^2 + reflect_dy^2 + reflect_dz^2);
          
    t = (z0 - Mirror_point(3))/unit_reflect_dz;
    arrive_point(1) = Mirror_point(1) + unit_reflect_dx * t;
    arrive_point(2) = Mirror_point(2) + unit_reflect_dy * t;
    arrive_point(3) = z0;
    
end