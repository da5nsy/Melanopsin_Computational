function melcomp_optimality(ss,prog,range)


load T_melanopsin T_melanopsin S_melanopsin; 
S_mel = SToWls(S_melanopsin);

[~,mloc] = max(T_melanopsin); %max location
Ip = S_mel(mloc); %Mel peak

if ~exist('range','var')
    range = [-50,5,50];
    disp(range);
end

for i=range(1):range(2):range(3)
    pc(i-range(1)+1) = melcomp(1,1,1,9,NaN,i);    
    ex(i-range(1)+1,:)=pc(i-range(1)+1).explained;
    if prog %progress
        disp(i)
    end
end

%plot(Ip+range(1):range(2):Ip+range(3),ex(1:range(2):end,1),'r','LineWidth',3)
%plot(Ip+range(1):range(2):Ip+range(3),ex(1:range(2):end,2),'g','LineWidth',3)
plot(Ip+range(1):range(2):Ip+range(3),ex(1:range(2):end,3),'b','LineWidth',3,'DisplayName','3rd PC explanation weight')
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
    
    legend
end

axis tight

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