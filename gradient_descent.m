function [polar_point,azi_point,P_max] = gradient_descent(P_NLoS,diff_polar,diff_azi,...
    A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr)
    
    l = 10^-1; % 梯度下降步长
    
    err = inf; % 定义初始误差
    
    Total_alpha = []; % 存储每次梯度下降的alpha值
    
    Total_beta = []; % 存储每次梯度下降的beta值
    
    Total_P = []; % 存储每次梯度下降的P_NLoS值
    
    i = 0;
    
    while i <= 1000
    
        i = i + 1;
        
        c_P = P_NLoS(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr) % 获取进行梯度变化前的P值
        diff_polar(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr)
        diff_azi(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr)
        alpha_RIS = alpha_RIS - l * diff_polar(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr); 
        
        beta_RIS = beta_RIS - l * diff_azi(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr);
        
        n_P = P_NLoS(A_PD,T,alpha_PD,alpha_RIS,alpha_source,beta_PD,beta_RIS,beta_source,f,hd,phi_semi,rho_mi,xd,xi_fov,xr,yd,yr,ys,zr)% 获取进行梯度变化后的P值
        
        err = abs(c_P-n_P); % 求误差
        
        Total_alpha = [Total_alpha;alpha_RIS];
        
        Total_beta = [Total_beta;beta_RIS];
        
        Total_P = [Total_P;n_P ];
        
        if err <= 10^-10 
        
            break;
        
        end
    
    end
    polar_point = alpha_RIS;
    azi_point = beta_RIS;
    P_max = n_P;
    
end

