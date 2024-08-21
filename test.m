% clear 
% close all
% figure
xlim([1,7])
ylim([1,7])
hold on
drawarrow([6,5],[6.5,5],'red','-','arrow',1.3)
text(6.52,5,'X','FontSize',14)
drawarrow([6,5],[6,5.8],'red','-','arrow',1.3)
text(5.9,5.9,'Y','FontSize',14)
drawarrow([6,5],[5.7,4.7],'red','-','arrow',1.3)
text(5.6,4.7,'Z','FontSize',14)
drawarrow([6,5],[6+0.6*cos(angle_PD(2)),5+0.25*sin(angle_PD(2))],'blue','-','arrow',1.3)
text(6.01+0.2*cos(angle_PD(2))+0.005,5+0.45*sin(angle_PD(2)-0.2),'n_{PD}=[1,-1,-1.41]','FontSize',14)
drawarrow([6,5],[6,4.55],'black',':','line',2)
drawarrow([6,5],[6.4,5.4],'black',':','line',2)
drawarrow([6+0.6*cos(angle_PD(2)),5+0.25*sin(angle_PD(2))],[6,5+0.25*sin(angle_PD(2))],'black',':','line',2)
text(5.88,5+0.25*sin(angle_PD(2)),'-1','FontSize',14)
drawarrow([6+0.6*cos(angle_PD(2)),5+0.25*sin(angle_PD(2))],[6.27,5.25],'black',':','line',2)
text(6.1,5.4,'-1.41','FontSize',14)
plot(coor_led(1),coor_led(2),'o','MarkerSize',6, 'MarkerFaceColor','red','LineWidth',1);
text(coor_led(1),coor_led(2)+0.2,'LED','FontSize',14)
xlabel('房间地面的X坐标');
ylabel('房间地面的Y坐标');
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'FontSize',16);