function melcomp_optimality(ss,prog,t_range)

try
    nargin;
catch
    %If we're not inside a function, clear everything and make a full-screen figure
    clear, clc, close all
    figure('units','normalized','outerposition',[0 0 1 1]), hold on
    
    %default values
    ss=1;
    prog=1;
    t_range = [-200,5,400];
    disp('Using default values');
end

%Compute melanopic peak sensitivity
load T_melanopsin T_melanopsin S_melanopsin; 
S_mel = SToWls(S_melanopsin);
[~,mloc] = max(T_melanopsin); %max location
Ip = S_mel(mloc); %Mel peak

t_range_ex=t_range(1):t_range(2):t_range(3);

for i=t_range_ex
    pc(i-t_range(1)+1) = melcomp_2(1,1,1,9,NaN,i); %principal components    
    ex(i-t_range(1)+1,:) = pc(i-t_range(1)+1).explained; %explained
    if prog %progress
        disp(i)
    end
end
ex_c = ex(1:t_range(2):end,3); %ex_condensed
[~,locs] = findpeaks(ex_c);
X_ax = Ip+t_range(1):t_range(2):Ip+t_range(3);
AoI = X_ax(locs); %Areas of interest

%Plot principal components
plt_firstAndSecondPC = 0;
plot(Ip+t_range_ex,ex_c,'b','LineWidth',3,'DisplayName','3rd PC weight')
if plt_firstAndSecondPC
    plot(Ip+t_range_ex,ex(1:t_range(2):end,1),'r','LineWidth',3,'DisplayName','1st PC weight')
    plot(Ip+t_range_ex,ex(1:t_range(2):end,2),'g','LineWidth',3,'DisplayName','2nd PC weight')
end

plt_melPeakLine = 0;
if plt_melPeakLine
    y=ylim; %Seems to jump around for some reason if I don't lock it down
    plot([Ip,Ip],[min(y),max(y)],'k--','DisplayName','Nominal melanopic peak sensitivity')
end

plt_thirdPCPeak = 0;
if plt_thirdPCPeak
    [~,m3] = max(ex(:,3)); %max third principal component
    plot([m3+Ip+range(1)-1,m3+Ip+range(1)-1],[min(y),max(y)],'DisplayName','3rd PC max')
end

if ss %spectral sensitvities
    load T_cones_sp T_cones_sp S_cones_sp
    load T_rods T_rods S_rods
    plot(SToWls(S_cones_sp),T_cones_sp(3,:),'b:','DisplayName','s-cone')
    plot(SToWls(S_melanopsin),T_melanopsin,':','DisplayName','melanopsin')
    plot(SToWls(S_rods),T_rods,':','DisplayName','rhodopsin')
    plot(SToWls(S_cones_sp),T_cones_sp(2,:),'g:','DisplayName','m-cone')
    plot(SToWls(S_cones_sp),T_cones_sp(1,:),'r:','DisplayName','l-cone')    
end

plt_peakLines = 0;
if plt_peakLines
    legend('AutoUpdate','off')
    for i = 1:length(locs)
        plot([AoI(i),AoI(i)],[min(ylim),pks(i)],'k')
    end
    % figure, plot(X_ax,ex_c) %safety check
end

legend('Location','best')
axis tight

%% See what the data looks like in 3 dimensions at the peaks pulled out by the previous calcs

plt_viz = 0;

if plt_viz
    for i = 1:length(locs)
        figure,
        melcomp_2(1,1,1,9,'3D',t_range_ex(locs(i)));
        title(Ip+t_range_ex(locs(i)))
    end
end

%%
% %Old way which plots the changing shape 
% for i=range(1):range(2)
%     figure('units','normalized','outerposition',[0 0 1 1]), hold on
%     subplot(1,2,1)
%     melcomp_2(1,1,1,9,'3D',i);
%     view(0,0)
%     zlim([0 1])
%     
%     subplot(1,2,2)
%     melcomp_2(1,1,1,9,'3D',i);
%     view(90,0)
%     zlim([0 1])
%     
%     drawnow
%     if i~=range(2)
%         clf
%     end
% end

end