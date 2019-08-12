function [sf_l,sf_s] = calcsf(lsri, l_cal_range, s_cal_range, plt, wholeset)

% This function calculates optimal scaling factors for k
% Modified from melcomp_6_calcsf.m

switch nargin
    case 0
        error('This function needs lsri at minimum')
    case 1
        l_cal_range = 0.35:0.0001:0.7;
        s_cal_range = -0.35:0.0001:-0.15;
        plt.iterations = 0;
        plt.sfs = 0;
        wholeset = 0;
end

if isfield(plt,'iterations') && plt.iterations % check it exists AND is turned on
    figure, hold on
end


%%

cal_val_std = zeros(2,size(l_cal_range,2)); %assuming the two ranges are the same length
for i=1:2 % for both l and s
    if i==1
        range = l_cal_range;
    else
        range = s_cal_range;
    end
    for cal_val = range
        if isfield(plt,'iterations') && plt.iterations
            cla reset
            if i == 1
                scatter(lsri(1,:)+cal_val*lsri(4,:),lsri(2,:),...
                    'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6) % would be nice to add colour
            else
                scatter(lsri(1,:),lsri(2,:)+cal_val*lsri(4,:),...
                    'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
            end
            title(cal_val)
            drawnow
            %pause(0.1)
        end
        if ~wholeset
            if i==1
                cal_val_std(i,find(range == cal_val)) = mean(std(squeeze(lsri(1,:,:)+cal_val*lsri(4,:,:))'));
            else
                cal_val_std(i,find(range == cal_val)) = mean(std(squeeze(lsri(2,:,:)+cal_val*lsri(4,:,:))'));
            end
        else
            if i==1
                cal_val_std(i,find(range == cal_val)) = std(lsri(1,:)+cal_val*lsri(4,:));
            else
                cal_val_std(i,find(range == cal_val)) = std(lsri(2,:)+cal_val*lsri(4,:));
            end
        end
    end
end

[~,cal_val_std_minloc] = min(cal_val_std,[],2);

sf_l = l_cal_range(cal_val_std_minloc(1));
sf_s = s_cal_range(cal_val_std_minloc(2));

%%

d = diff(cal_val_std);
if or(~(and(d(1)<0,d(end)>0)),(isfield(plt,'sfs') && plt.sfs)) % If you've requested this plot, or if there's an error...
    figure, 
    plot(l_cal_range,cal_val_std)
    xlabel('k')
    ylabel('SD')
    legend('k1','k2')
    if ~(and(d(1)<0,d(end)>0))
        error('Minima for l not reached. Modify search range')
    end
end

end

