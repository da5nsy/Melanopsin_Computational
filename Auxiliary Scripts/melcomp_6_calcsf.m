function [sf_l,sf_s] = melcomp_6_calcsf(lsri, l_cal_range, s_cal_range,plot_iterations,pltc_alt)

switch nargin
    case 0
        error('This function needs lsri at minimum')
    case 1
        l_cal_range = 0.35:0.0001:0.7;
        s_cal_range = -0.35:0.0001:-0.15;
        plot_iterations = 0;
    case 3
        plot_iterations = 0;
end

if plot_iterations
    figure, hold on
end

%% 

l_cal_val_std = [];
for l_cal_val = l_cal_range
    if plot_iterations
        cla reset
        scatter(lsri(1,:)+l_cal_val*lsri(4,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
        title(l_cal_val)
        drawnow
        pause(0.1)
    end
    l_cal_val_std = [l_cal_val_std, mean(std(squeeze(lsri(1,:,:)+l_cal_val*lsri(4,:,:))'))];
end

[~,l_cal_val_std_minloc] = min(l_cal_val_std);
sf_l = l_cal_range(l_cal_val_std_minloc);

%%

s_cal_val_std = [];
for s_cal_val = s_cal_range
    if plot_iterations
        cla reset
        scatter(lsri(1,:),lsri(2,:)+s_cal_val*lsri(4,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
        title(s_cal_val)
        drawnow
        pause(0.1)
    end
    s_cal_val_std = [s_cal_val_std, mean(std(squeeze(lsri(2,:,:)+s_cal_val*lsri(4,:,:))'))];
end


[~,s_cal_val_std_minloc] = min(s_cal_val_std);
sf_s = s_cal_range(s_cal_val_std_minloc);

%%

d = diff(l_cal_val_std);
if ~(and(d(1)<0,d(end)>0))
    figure, plot(l_cal_range,l_cal_val_std)
    error('Minima for l not reached. Modify search range')    
end

d = diff(s_cal_val_std);
if ~(and(d(1)<0,d(end)>0))
    figure, plot(s_cal_range,s_cal_val_std)
    error('Minima for s not reached. Modify search range')
end

%% Force show the above figures, for debugging

% figure, plot(l_cal_range,l_cal_val_std)
% figure, plot(s_cal_range,s_cal_val_std)

end

