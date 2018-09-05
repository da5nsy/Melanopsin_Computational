function melcomp_optimality(ss,prog,t_range)

try
    nargin;
catch
    %If we're not inside a function, clear everything and make a full-screen figure
    clear, clc, close all
    figure('units','normalized','outerposition',[0 0 1 1]), hold on
    
    %default values
    ss=1;
    prog=0;
    t_range = [-200,10,400];
    disp('Using default values');
end

load T_melanopsin T_melanopsin S_melanopsin; 
S_mel = SToWls(S_melanopsin);

[~,mloc] = max(T_melanopsin); %max location
Ip = S_mel(mloc); %Mel peak

for i=t_range(1):t_range(2):t_range(3)
    pc(i-t_range(1)+1) = melcomp(1,1,1,9,NaN,i);    
    ex(i-t_range(1)+1,:)=pc(i-t_range(1)+1).explained;
    if prog %progress
        disp(i)
    end
end

%plot(Ip+range(1):range(2):Ip+range(3),ex(1:range(2):end,1),'r','LineWidth',3)
%plot(Ip+range(1):range(2):Ip+range(3),ex(1:range(2):end,2),'g','LineWidth',3)
plot(Ip+t_range(1):t_range(2):Ip+t_range(3),ex(1:t_range(2):end,3),'b','LineWidth',3,'DisplayName','3rd PC explanation weight')
y=ylim; %Seems to jump around for some reason if I don't lock it down
plot([Ip,Ip],[min(y),max(y)],'k--','DisplayName','Nominal melanopic peak sensitivity')

% [~,m3] = max(ex(:,3)); %max third principal component
% plot([m3+Ip+range(1)-1,m3+Ip+range(1)-1],[min(y),max(y)],'DisplayName','3rd PC max')

if ss %spectral sensitvities
    load T_cones_sp
    load T_rods
    plot(SToWls(S_cones_sp),T_cones_sp(3,:),'b:','DisplayName','s-cone')
    plot(SToWls(S_melanopsin),T_melanopsin,':','DisplayName','melanopsin')
    plot(SToWls(S_rods),T_rods,':','DisplayName','rhodopsin')
    plot(SToWls(S_cones_sp),T_cones_sp(2,:),'g:','DisplayName','m-cone')
    plot(SToWls(S_cones_sp),T_cones_sp(1,:),'r:','DisplayName','l-cone')
    
    legend('Location','best')
end

axis tight

%%
ex_c = ex(1:t_range(2):end,3); %ex_condensed
[pks,locs] = findpeaks(ex_c);

%!!!!!!!!!!
%locs in reference to original index (aka AoI) = ??;
%%
viz = 1;
 
if viz
    AoI = []; %Areas of interest
    for i = AoI
    melcomp(1,1,1,9,'3D',i);
    end
end

%%
% %Old way which plots the changing shape 
% for i=range(1):range(2)
%     figure('units','normalized','outerposition',[0 0 1 1]), hold on
%     subplot(1,2,1)
%     melcomp(1,1,1,9,'3D',i);
%     view(0,0)
%     zlim([0 1])
%     
%     subplot(1,2,2)
%     melcomp(1,1,1,9,'3D',i);
%     view(90,0)
%     zlim([0 1])
%     
%     drawnow
%     if i~=range(2)
%         clf
%     end
% end

end