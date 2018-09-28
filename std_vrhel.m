load sur_vrhel.mat

figure, hold on
plot(SToWls(S_vrhel),std(sur_vrhel'))
plot(SToWls(S_vrhel),std(sur_vrhel')./mean(sur_vrhel'),'g')

legend

load T_cones_sp.mat T_cones_sp S_cones_sp
plot(SToWls(S_cones_sp),T_cones_sp)

%%

clc, clear

load sur_vrhel.mat
load B_cieday.mat

D65 = GenerateCIEDay(6500,B_cieday);
D65_r = SplineSpd(S_cieday,D65,S_vrhel);

for i=1:size(sur_vrhel,2)
    sur_vrhel(:,i) = sur_vrhel(:,i).*D65_r;
end

%figure, hold on
%plot(SToWls(S_vrhel),std(sur_vrhel'))
plot(SToWls(S_vrhel),std(sur_vrhel')./mean(sur_vrhel'),'r--')

