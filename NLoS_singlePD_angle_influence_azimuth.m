clear
close all

pd_polar_save = pi/2:pi/8:pi;
for i = 1:length(pd_polar_save)
    NLoS_singlePD_plot_azimuth(pd_polar_save(i));
end