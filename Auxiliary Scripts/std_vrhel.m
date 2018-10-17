% Where would you put (spectrally) a sensor to best judge how much light
% there is on an object?
% I guess you'd want somewhere where there was minimal intrusion into
% variability from differential reflectances, aka where variability in
% reflectance was lowest vs power.
% (There might be an additional constraint that you'd want enough power to
% avoid SNR issues)

% TO DO LIST
% ----------

% Add back in D65 from previous version. You'd want to sample the area of
% the spectrum where there was the most power, also.

% Try for Foster images

load sur_vrhel.mat

figure, hold on
plot(SToWls(S_vrhel),std(sur_vrhel'),'DisplayName','SD')
plot(SToWls(S_vrhel),mean(sur_vrhel'),'DisplayName','Mean')
plot(SToWls(S_vrhel),std(sur_vrhel')./mean(sur_vrhel'),'DisplayName','SD/Mean')
% plot(D65) !!!!!!!!!!

legend('Autoupdate','off')

load T_cones_sp.mat T_cones_sp S_cones_sp
plot(SToWls(S_cones_sp),T_cones_sp,':')
