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
    t_range = [-100,5,250];
    disp('Using default values');
end

%Compute melanopic peak sensitivity
load T_melanopsin T_melanopsin S_melanopsin; 
S_mel = SToWls(S_melanopsin);
[~,mloc] = max(T_melanopsin); %max location
Ip = S_mel(mloc); %Mel peak

t_range_ex=t_range(1):t_range(2):t_range(3);

for i=t_range_ex
    pc(i-t_range(1)+1) = melcomp_2('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','Z_ax',9,'mel_offset',i); %principal components    
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
plt_firstAndSecondPC = 1;
if plt_firstAndSecondPC
    plot(Ip+t_range_ex,ex(1:t_range(2):end,1)/max(ex(1:t_range(2):end,1)),...
        'r','LineWidth',3,'DisplayName',['1st PC weight /' num2str(max(ex(1:t_range(2):end,1)))])
    plot(Ip+t_range_ex,ex(1:t_range(2):end,2)/max(ex(1:t_range(2):end,2)),...
        'g','LineWidth',3,'DisplayName',['2nd PC weight /' num2str(max(ex(1:t_range(2):end,2)))])
end
plot(Ip+t_range_ex,ex(1:t_range(2):end,3)/max(ex(1:t_range(2):end,3)),...
    'b','LineWidth',3,'DisplayName',['3rd PC weight /' num2str(max(ex(1:t_range(2):end,3)))])

plt_melPeakLine = 1;
if plt_melPeakLine
    y=ylim; %Seems to jump around for some reason if I don't lock it down
    plot([Ip,Ip],[min(y),max(y)],'k--','DisplayName','Nominal melanopic peak sensitivity')
end

plt_thirdPCPeak = 1;
if plt_thirdPCPeak
    [~,m3] = max(ex(:,3)); %max third principal component
    plot([m3+Ip+t_range(1)-1,m3+Ip+t_range(1)-1],[min(y),max(y)],'DisplayName','3rd PC max')
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

% %broken, and can't work out what it is meant to do
% plt_peakLines = 0;
% if plt_peakLines
%     legend('AutoUpdate','off')
%     for i = 1:length(locs)
%         plot([AoI(i),AoI(i)],[min(ylim),pks(i)],'k')
%     end
%     % figure, plot(X_ax,ex_c) %safety check
% end

legend('Location','best')
axis tight

% %% See what the data looks like in 3 dimensions at the peaks pulled out by the previous calcs
% 
% plt_viz = 1;
% 
% if plt_viz
%     for i = 1:length(locs)       
%         melcomp_2('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','Z_ax',9,'plt','3D','mel_offset',t_range_ex(locs(i)));
%         title(Ip+t_range_ex(locs(i)))
%     end
% end
% 
% %%
% %Old way which plots the changing shape 
% %need to turn 'figure' command in melcomp_2 plt_3D off to get this to
% %display as planned
% 
% t_range = [-100,2,200];
% 
% figure('units','normalized','outerposition',[0.3 0.3 0.5 0.5]), hold on
% 
% for i=t_range(1):t_range(2):t_range(3)
%     subplot(1,2,1)
%     melcomp_2('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','Z_ax',9,'plt','3D','mel_offset',i);
%     view(20,10)
%     camproj('perspective')
%     zlim([0 3])
%     text(0.1,0.1,num2str(i))
%     rmv_lbls
%     
%     subplot(1,2,2)
%     melcomp_2('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','Z_ax',9,'plt','3D','mel_offset',i);
%     view(22,10)
%     camproj('perspective')
%     zlim([0 3])
%     text(0.1,0.1,num2str(i))
%     rmv_lbls
%    
%     drawnow
%     if i == t_range(1)
%         waitforbuttonpress
%     end
%     if i ~= t_range(3)
%         clf
%     end
%     
% end
% 
% %%
% 
% t_range = [-100,20,200];
% lt = length(t_range(1):t_range(2):t_range(3));
% counter = 1;
% 
% figure('units','normalized','outerposition',[0.3 0.3 0.5 0.5]), hold on
% 
% for i=t_range(1):t_range(2):t_range(3)
%     subplot(ceil(sqrt(lt)),ceil(sqrt(lt)),counter)
%     melcomp_2('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','Z_ax',9,'plt','3D','mel_offset',i);
%     view(-170,45)
%     xlim([0.5 1])
%     zlim([0 1])
%     text(0.6,0.1,num2str(i))
%     rmv_lbls
%    
%     drawnow  
%     counter = counter +1;
%     
% end

end