%clear, close all

figure,
melcomp_2('SPD','Granada_sub','SRF','Vrhel_nat_extended','plt','3D');
axis equal

xlim([-3 3])
ylim([-3 3])
zlim([-3 3])

xticks([-3 0 3])
yticks([-3 0 3])
zticks([-3 0 3])

grid on

xlabel('l*')
ylabel('s*')
zlabel('i*')

view(2)