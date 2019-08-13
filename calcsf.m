%function [sf_l,sf_s,minSD_l,minSD_s] = calcsf(lsri, l_cal_range, s_cal_range, plt, wholeset, map)
function [sf_l,sf_s,minSD_l,minSD_s] = calcsf(lsri,varargin)

% This function calculates optimal scaling factors for k
% Modified from melcomp_6_calcsf.m

% Default behavoir assumes you want to minimise SD for individual
% groupings, turn on 'wholeset' flag if you want to minimise variation
% across the whole set.

%% Parse inputs

defaultplt.iterations = 0;
defaultplt.sfs = 0;
defaultmap = repmat(1:size(lsri,2),size(lsri,3),1)';

p = inputParser;

addParameter(p,'l_cal_range',linspace(-5,5,401));
addParameter(p,'s_cal_range',linspace(-5,5,401));
addParameter(p,'plt',defaultplt);
addParameter(p,'wholeset',0);
addParameter(p,'map',defaultmap);

parse(p,varargin{:});

l_cal_range = p.Results.l_cal_range;
s_cal_range = p.Results.s_cal_range;
plt = p.Results.plt;
wholeset = p.Results.wholeset;
map = p.Results.map;

if ~isfield(plt,'iterations')
    plt.iterations = defaultplt.iterations;
end

if ~isfield(plt,'sfs')
    plt.sfs = defaultplt.sfs;
end

if length(l_cal_range) ~= length(s_cal_range)
    error('This function needs l_cal_range and s_cal_range to be the same length. I know that this is not ideal. Sorry.')
end

%%
if plt.iterations
    figure, hold on
end

cal_val_std = zeros(2,size(l_cal_range,2)); %assuming the two ranges are the same length
for ls=1:2 % for both l and s
    if ls==1
        range = l_cal_range;
    else
        range = s_cal_range;
    end
    for cal_val_i = 1:length(range)
        if plt.iterations
            cla reset
            if ls == 1
                scatter(lsri(1,:)+range(cal_val_i)*lsri(4,:),lsri(2,:),...
                    'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6) % would be nice to add colour
            else
                scatter(lsri(1,:),lsri(2,:)+range(cal_val_i)*lsri(4,:),...
                    'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
            end
            title(range(cal_val_i))
            drawnow
            %pause(0.1)
            
        end
        if ~wholeset
            for surfI = 1:size(lsri,2)%surfindex
                SDstore(surfI) = std(lsri(ls,map == surfI)+range(cal_val_i)*lsri(4,map == surfI));
            end
            cal_val_std(ls,cal_val_i) = mean(SDstore);
        else
            cal_val_std(ls,cal_val_i) = std(lsri(ls,:)+range(cal_val_i)*lsri(4,:));
        end
    end
end

[minSD,cal_val_std_minloc] = min(cal_val_std,[],2);

minSD_l = minSD(1);
minSD_s = minSD(2);

sf_l = l_cal_range(cal_val_std_minloc(1));
sf_s = s_cal_range(cal_val_std_minloc(2));

%%

d1 = diff(cal_val_std(1,:));
d2 = diff(cal_val_std(2,:));
% haven't checked the below, use previous version if it doesn't work
if or(plt.sfs,(or(~(and(d1(1)<0,d1(end)>0)),~(and(d2(1)<0,d2(end)>0))))) % If you've requested this plot, or if there's an error...
    figure,
    plot(l_cal_range,cal_val_std)
    xlabel('k')
    ylabel('SD')
    legend('k1','k2')
    if or(~(and(d1(1)<0,d1(end)>0)),~(and(d2(1)<0,d2(end)>0)))
        error('Minima for l not reached. Modify search range')
    end
end

end

