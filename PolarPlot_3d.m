function [] = PolarPlot_3d(figure_name,polar,azi,distance,coor_PD,x_name,y_name,z_name,title_name)

[polar_3d, azi_3d] = meshgrid(unique(polar), unique(azi));
distance_3d = reshape(distance,length(polar),length(azi));
% figure
% trisurf(delaunay(polar_3d(:), azi_3d(:)), polar_3d(:), azi_3d(:), distance_3d(:), 'FaceColor', 'interp', 'EdgeColor', 'k');
x = distance_3d.*sin(polar_3d).*cos(azi_3d)+coor_PD(1);
y = distance_3d.*sin(polar_3d).*sin(azi_3d)+coor_PD(2);
z = -distance_3d.*cos(polar_3d);

figure(figure_name)
plot3(x,y,z);
% xticks([0 pi/8 2*pi/8 3*pi/8 4*pi/8 5*pi/8 6*pi/8 7*pi/8 8*pi/8])
% xticklabels({'0','\frac{\pi}{8}','\frac{\pi}{4}','\frac{3\pi}{8}','\frac{\pi}{2}','\frac{5\pi}{8}','\frac{3\pi}{4}','\frac{7\pi}{8}','\pi'})
% xticks([0 pi/4 2*pi/4 3*pi/4 4*pi/4 5*pi/4 6*pi/4 7*pi/4 8*pi/4])
% xticklabels({'0','\frac{\pi}{4}','\frac{\pi}{2}','\frac{3\pi}{4}','\pi','\frac{5\pi}{4}','\frac{3\pi}{2}','\frac{7\pi}{4}','2\pi'})
% 
xlabel(x_name)
ylabel(y_name)
zlabel(z_name)
title(title_name)

end