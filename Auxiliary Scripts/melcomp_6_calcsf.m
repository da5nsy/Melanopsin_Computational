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
    end
    l_cal_val_std = [l_cal_val_std, std(lsri(1,:)+l_cal_val*lsri(4,:))];
end

[~,l_cal_val_std_minloc] = min(l_cal_val_std);
sf_l = l_cal_range(l_cal_val_std_minloc);

%%

s_cal_val_std = [];
for s_cal_val = s_cal_range
    if plot_iterations
        cla reset
        scatter(lsri(1,:),lsri(2,:)+s_cal_val*lsri(4,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
        title(l_cal_val)
        drawnow
    end
    s_cal_val_std = [s_cal_val_std, std(lsri(2,:)+s_cal_val*lsri(4,:))];
end


[~,s_cal_val_std_minloc] = min(s_cal_val_std);
sf_s = s_cal_range(s_cal_val_std_minloc);

%%
% figure, plot(l_cal_range,l_cal_val_std)
% figure, plot(s_cal_range,s_cal_val_std)

end

