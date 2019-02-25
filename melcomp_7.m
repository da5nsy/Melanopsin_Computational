

%%

load T_cones_ss2.mat

%%
plot3(T_cones_ss2(1,:),T_cones_ss2(2,:),T_cones_ss2(3,:),'k.')

xlabel('L')
ylabel('M')
zlabel('S')

axis equal

xlim([0,1])
ylim([0,1])
zlim([0,1])

cleanTicks
grid on