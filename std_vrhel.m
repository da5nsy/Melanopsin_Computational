load sur_vrhel.mat

figure, hold on
plot(SToWls(S_vrhel),std(sur_vrhel'))
plot(SToWls(S_vrhel),std(sur_vrhel')./mean(sur_vrhel'),'g')

legend

load T_cones_sp.mat T_cones_sp S_cones_sp
plot(SToWls(S_cones_sp),T_cones_sp)
