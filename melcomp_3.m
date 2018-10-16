try
    nargin;
catch
    clear, clc, close all
end

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Project Overview Presentation';

ff = '-dtiff'; %file format
p  = 1; %print? (aka save?), set to 1 to begin saving

plot_where = [800,50];
plot_size  = [800,375];

%% Start the figure

figure('Position',[plot_where plot_size],'defaultLineLineWidth',2)
hold on
set(gca, 'FontSize', 16)

%% SPDs (Spectral Power Distributions)

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

% try load SPD.mat addI %bit slow so I saved it out
% catch
%     % Additional SPD info (source: personal correspondance with J. Hernández-Andrés)
%     [addI.NUM,addI.TXT,addI.RAW] = xlsread('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\add_info.xlsx');
%     for i=1:length(T_SPD) %a lot of stupid code just to get the date and time out
%         addI.t(i) = datetime(...
%             [char(addI.RAW(2,2)),' ',char(days(cell2mat(addI.RAW(2,3))),'hh:mm')],...
%             'InputFormat','dd/MM/uuuu HH:mm');
%     end
%     addI.el = addI.NUM(:,4); %elevation
%     addI.az = addI.NUM(:,5); %azimuth
% end

% T_SPD = T_SPD(:,addI.el>30);
% addI.el = addI.el(addI.el>30);

plt_SPD = 1;
if plt_SPD
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([0 2]); yticks([0,1,2]);
    xlabel('Wavelength'); ylabel('Power');
    
    plot(SToWls([380,5,81]),T_SPD(17:97,103)) %single spd
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls([380,5,81]),T_SPD(17:97,:)) %all spd
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%% SRFs (Spectral Reflectance Functions)

load sur_vrhel
T_SRF = sur_vrhel(:,[1:44,65,69,81:154]); %natural only
S_SRF = S_vrhel;

plt_refs = 1;
if plt_refs
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([0 1]); yticks([0,1]);
    xlabel('Wavelength'); ylabel('Reflectance');
    
    plot(SToWls(S_SRF),T_SRF(:,109)) %pear
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_SRF),T_SRF) %all
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%% SSF (Spectral Sensitivity Functions)

% change for ss at some point in future (?)
load T_cones_sp T_cones_sp S_cones_sp
T_SSF = T_cones_sp;
S_SSF = S_cones_sp;
clear T_cones_sp S_cones_sp

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin, T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0];
clear S_melanopsin T_melanopsin

plt_spec_sens = 1;
if plt_spec_sens
    figure(1)
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([0 1]); yticks([0,1]);
    ylabel('Sensitivity')
    plot(SToWls(S_SSF),T_SSF)
    plot(SToWls(S_rods),T_rods)
    plot(SToWls(S_mel),T_mel)
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%% Interpolate functions to match range and interval

%reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1)]),max([S_SPD(2),S_SRF(2),S_SSF(2)]),min([S_SPD(3),S_SRF(3),S_SSF(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is
%S_sh = [400,5,61]; % brute force way to do something like the variableweights thing

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_SRF = SplineSrf(S_SRF,T_SRF,S_sh,1);
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)'; %extend? !!!!!!!!!!
T_mel  = SplineCmf(S_mel,T_mel,S_sh)'; %extend? !!!!!!!!!!
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

% compute peaks of spectral sensitivities
[~,ms_ind] = max(T_LMSRI);
S_full = SToWls(S_sh);
ms_val = S_full(ms_ind);

%% PCA

% PCA of spectral Power distributions
% compute pca variable weight
vw = 1;
if vw == 0 %no variable weights
    pc_p.variableweights = ones(81,1);
elseif vw == 1 %S + L and max in between
    pc_p.variableweights = T_SSF(:,1)+T_SSF(:,3);
    [~,t_p_locs] = findpeaks(pc_p.variableweights);
    pc_p.variableweights(t_p_locs(1):t_p_locs(2)) = max(pc_p.variableweights);
    %plot(pc.variableweights)
elseif vw == 2 %S+L+I
    pc_p.variableweights = T_SSF(:,1)+T_SSF(:,3)+T_mel(:,1);
end

[pc_p.coeff, pc_p.score, pc_p.latent, pc_p.tsquared, pc_p.explained, pc_p.mu] = pca(T_SPD','VariableWeights',pc_p.variableweights);

plt_pc_p = 1;
if plt_pc_p
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([-0.55 0.55]); yticks([-0.5,0,0.5]);
    ylabel('Power')
    %legend({'PC1','PC2','PC3'},'Location','Southwest')
    
    % highlights negative space
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    
    %normalised PCs:
    %plot(SToWls(S_SPD),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))
    
    plot(SToWls(S_SPD),pc_p.coeff(:,1))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_SPD),pc_p.coeff(:,2))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_SPD),pc_p.coeff(:,3))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    hold on
    ylim([0 2]); yticks([0,1,2]);
    
    %reconstructed data
    recon_p(:,1) = pc_p.score(103,:)   * pc_p.coeff'        + pc_p.mu; % Reconstructed with all PCs
    recon_p(:,2) = pc_p.score(103,1:3) * pc_p.coeff(:,1:3)' + pc_p.mu; % Reconstructed with 3 PCs
    
    %plot(SToWls(S_SPD), T_SPD(:,103), 'b') % Original data
    %plot(SToWls(S_SPD), recon_p(:,1), 'g') % All PC
    plot(SToWls(S_SPD), recon_p(:,2), 'r') % 3 PC
    % legend({'Original Data','Reconstructed with 1:3'})
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    % figure, plot((pc_p.score*pc_p.coeff' + repmat(pc_p.mu,length(pc_p.score),1))') %all the data reconstructed
    
    % disp(pc_p.score(103,1:3)) %values for demo ref
end

[pc_r.coeff, pc_r.score, pc_r.latent, pc_r.tsquared, pc_r.explained, pc_r.mu] = pca(T_SRF','VariableWeights',pc_p.variableweights);
%inhereted variable weights from pc_p

plt_pc_r = 1;
if plt_pc_r
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([-0.55 0.55]); yticks([-0.5,0,0.5]);
    ylabel('Power')
    %legend({'PC1','PC2','PC3'},'Location','Southwest')
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    
    %plot(SToWls(S_SRF),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))
    
    plot(SToWls(S_SRF),pc_r.coeff(:,1))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_SRF),pc_r.coeff(:,2))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_SRF),pc_r.coeff(:,3))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    hold on
    ylim([0 1]); yticks([0,1]);
    
    recon_r(:,1) = pc_r.score(109,:)   * pc_r.coeff'        + pc_r.mu; % Reconstructed with all PCs
    recon_r(:,2) = pc_r.score(109,1:3) * pc_r.coeff(:,1:3)' + pc_r.mu; % Reconstructed with 3 PCs
    
    %plot(SToWls(S_SRF), T_SRF(:,109), 'b') % Original data
    %plot(SToWls(S_SRF), recon_r(:,1), 'g') % All PC
    plot(SToWls(S_SRF), recon_r(:,2), 'r') % 3 PC
    % legend({'Original Data','Reconstructed with 1:3'})
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    % figure, plot((pc_r.score*pc_r.coeff' + repmat(pc_r.mu,length(pc_r.score),1))') %all the data reconstructed
    
    % disp(pc_r.score(109,1:3)) %values for demo ref
end

% Show the effect of changing the weight of PC2
plt_pc2_demo = 1;
if plt_pc2_demo
    figure(1)
    cols = lines(3);
    xlim([380 780]); xticks(380:100:780);
    ylim([-0.55 0.55]); yticks([-0.5,0,0.5]);
    ylabel('Power');
    
    cla
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    
    plot(SToWls(S_SPD),pc_p.coeff(:,2),'Color',cols(2,:))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot([ms_val(3),ms_val(3)],[min(ylim),max(ylim)],'Color','k')
    plot([ms_val(5),ms_val(5)],[min(ylim),max(ylim)],'Color','k')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    plot(SToWls(S_SPD),pc_p.coeff(:,2)*-1,'Color',cols(2,:)*((-1+2)/4))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot([ms_val(3),ms_val(3)],[min(ylim),max(ylim)],'Color','k')
    plot([ms_val(5),ms_val(5)],[min(ylim),max(ylim)],'Color','k')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    for i = -1.5:0.5:1.5
        plot(SToWls(S_SPD),pc_p.coeff(:,2)*i,'Color',cols(2,:)*((i+3)/4))
    end
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end


%% Calculate colour signals

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i)); %change this to be more direct !!!!!!!!!!!!!!
    lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); % double check this !!!!!!!!!!!!!!
    lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); % double check this !!!!!!!!!!!!!!
end

%% Calculate lines of fit between S and I values

for i=1:size(LMSRI,3)
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
end 

f_m = permute(repmat(fit(:,1),1,1,size(LMSRI,2)),[2,3,1]); % Mx + c
f_c = permute(repmat(fit(:,2),1,1,size(LMSRI,2)),[2,3,1]); % mx + C

% - % Checking

% figure('defaultLineLineWidth',2), hold on
% xlabel('log(S)'); ylabel('log(I)');
% %axis equal
% 
% for i = 1:10:2600
%     %scatter3(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),ones(120,1)*i,'k','filled')
%         
%     x_t = linspace(min(log10(LMSRI(3,:,i))),max(log10(LMSRI(3,:,i))));
%     y_t = x_t*fit(i,1) + fit(i,2);
%     
%     x_t2 = linspace(min(log10(LMSRI(3,:,i))),max(log10(LMSRI(3,:,i))));
%     y_t2 = x_t*f_m(1,1,i) + f_c(1,1,i);
%     
%     plot3(x_t,y_t,ones(100,1)*i,'r')
%     %plot3(x_t2,y_t2,ones(100,1)*i,'g:')%     
% end
% 
% view(-43,0)
% 
% figure,
% scatter3(f_m(1,1,:),f_c(1,1,:),pc_p.score(:,2))
% xlabel('f_m'); ylabel('f_c');



%% Calc correlation between

nPC = 3; %Number of principal components

cs = cat(1,LMSRI,...
    lsri,...
    LMSRI(1,:,:)+LMSRI(2,:,:),...
    (0.6373*LMSRI(1,:,:))+(0.3924*LMSRI(2,:,:)),...
    LMSRI(4,:,:)+ LMSRI(5,:,:),...
    LMSRI(1,:,:)./LMSRI(2,:,:),...
    LMSRI(1,:,:)./LMSRI(3,:,:),...
    LMSRI(1,:,:)./LMSRI(4,:,:),...
    LMSRI(1,:,:)./LMSRI(5,:,:),...
    LMSRI(3,:,:)./LMSRI(4,:,:),...
    LMSRI(3,:,:)./LMSRI(5,:,:),...
    (LMSRI(3,:,:)+LMSRI(5,:,:))./(LMSRI(1,:,:)+LMSRI(2,:,:)),...
    f_m,...
    f_c...
    ); %colour signals

plt_lbls{1}  = 'L'; %writing out this way so that there's a quick reference as to which value of Z_ax does what
plt_lbls{2}  = 'M';
plt_lbls{3}  = 'S';
plt_lbls{4}  = 'R';
plt_lbls{5}  = 'I';
plt_lbls{6}  = 'l';
plt_lbls{7}  = 's';
plt_lbls{8}  = 'r';
plt_lbls{9}  = 'i';
plt_lbls{10} = 'L+M';
plt_lbls{11} = '(0.6373*L)+(0.3924*M)';
plt_lbls{12} = 'r + i';
plt_lbls{13} = 'L/M';
plt_lbls{14} = 'L/S';
plt_lbls{15} = 'L/R';
plt_lbls{16} = 'L/I';
plt_lbls{17} = 'S/R';
plt_lbls{18} = 'S/I';
plt_lbls{19} = '(S+I)/(L+M)';
plt_lbls{20} = 'S/I fit m';
plt_lbls{21} = 'S/I fit c';

c = zeros(size(cs,1),nPC);
for j=1:nPC
    for i=1:size(cs,1)
        c(i,j) = corr(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j));
    end
end
c=abs(c);

c_norm = c;
for i=1:j %leftover, be careful
    c_norm(:,i) = c_norm(:,i) - min(c_norm(:,i));
    c_norm(:,i) = c_norm(:,i)/max(c_norm(:,i));
end

plt_viz = 1;
if plt_viz
    figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
    hold on
    
    subplot(1,2,1)
    imagesc(c)
    colorbar
    title('correlation between signal and PC weight')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:size(cs,1));
    set(gca, 'YTickLabel', plt_lbls);
    colormap('gray');
    set(gca, 'FontSize', 16)
    
    subplot(1,2,2)
    imagesc(c_norm)
    colorbar
    title('normalised')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:size(cs,1));
    set(gca, 'YTickLabel', plt_lbls);
    colormap('gray');   
    set(gca, 'FontSize', 16)
    
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

% plot scatters with backrounds the colours of the strength of correlation

plt_sca = 0;
if plt_sca
    figure,
    for i=1:size(cs,1)
        for j=1:nPC
            subplot(size(cs,1),nPC,j+((i-1)*nPC))
            scatter(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j),'r.')
            if j == 1
                ylabel(plt_lbls{i},'rotation',0)
            end
            set(gca,'YTickLabel',[])
            set(gca,'XTickLabel',[])
            set(gca,'Color',repmat(c_norm(i,j),3,1))
            axis tight
        end
    end
end


%% Why S and I when larger intervals are available?

% Correlation image across reflectance spectra
figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
hold on
set(gca, 'FontSize', 16)

sur_vrhel_n = sur_vrhel(:,[1:44,65,69,81:154]);
sur_vrhel_c = corr(sur_vrhel_n');
surf(sur_vrhel_c*100,'EdgeColor','none')
axis image
colormap gray
S_vrhel_f = SToWls(S_vrhel);
xticks(6:50:200)
yticks(6:50:200)
set(gca,'XTickLabel',S_vrhel_f(xticks))
set(gca,'YTickLabel',S_vrhel_f(yticks))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot3([26,26],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','b')
plot3([50,50],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','k')
plot3([min(xlim),max(xlim)],[26,26],[max(zlim),max(zlim)],'Color','b')
plot3([min(xlim),max(xlim)],[50,50],[max(zlim),max(zlim)],'Color','k')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

plot3([88,88],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','r')
plot3([78,78],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color','g')
plot3([58,58],[min(ylim),max(ylim)],[max(zlim),max(zlim)],'Color',[0.5,0.5,0.5])
plot3([min(xlim),max(xlim)],[88,88],[max(zlim),max(zlim)],'Color','r')
plot3([min(xlim),max(xlim)],[78,78],[max(zlim),max(zlim)],'Color','g')
plot3([min(xlim),max(xlim)],[58,58],[max(zlim),max(zlim)],'Color',[0.5,0.5,0.5])
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%-% Plot normalised spectral reflectances
figure(1)
cla
hold on
xlim([390 730]); xticks(400:100:700);
ylim([0 10]); yticks([0,10]);
xlabel('Wavelength'); ylabel('Normalised  Reflectance');

plot([440,440],[min(ylim),max(ylim)],'b:')
plot([488,488],[min(ylim),max(ylim)],'k:')

plot(SToWls(S_vrhel),sur_vrhel_n./sur_vrhel_n(26,:))
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% 

%L against PC1
figure(1)
cla
scatter(squeeze(mean(cs(1,:,:),2)),pc_p.score(:,1),'k.')
xlabel(plt_lbls{1})
ylabel('PC 1')
xlim('auto')
xticks('auto')
ylim('auto')
yticks('auto')

f_L1 = polyfit(squeeze(mean(cs(1,:,:),2)),pc_p.score(:,1),1);
x = linspace(min(squeeze(mean(cs(1,:,:),2))),max(squeeze(mean(cs(1,:,:),2))));
y = x*f_L1(1)+f_L1(2);
plot(x,y)

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

% S/I against PC2
cla
scatter(squeeze(mean(cs(18,:,:),2)),pc_p.score(:,2),'k.')
xlabel(plt_lbls{18})
ylabel('PC 2')

f_SI2 = polyfit(squeeze(mean(cs(18,:,:),2)),pc_p.score(:,2),1);
x = linspace(min(squeeze(mean(cs(18,:,:),2))),max(squeeze(mean(cs(18,:,:),2))));
y = x*f_SI2(1)+f_SI2(2);
plot(x,y)

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% 3D plot - basic (same as previous but without line and with a third dimension)
cla
scatter3(squeeze(mean(cs(18,:,:),2)),pc_p.score(:,2),pc_p.score(:,1),'k.')
%set(gca,'Color',repmat(c_norm(18,2),3,1))
xlabel(plt_lbls{18})
ylabel('PC 2')
zlabel('PC 1')

% how to visualize in 2D? !!!!!!!!!!!!!!!!!1
% if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% 3D plot - log(S/I) against (normalised(PC1)).^(1/4)

pc1n = (pc_p.score(:,1)-min(pc_p.score(:,1))).^(1/4);

cla
scatter3(log2(squeeze(mean(cs(18,:,:),2))),pc_p.score(:,2),pc1n,'b.')
xlabel('log2(S/I)')
ylabel('PC2')
zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')

% how to visualize in 2D? !!!!!!!!!!!!!!!!!1
% if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Plot showing segmentation by PC1

plt_seg = 1;

NoD = 30; %nnumber of divisions
m1 = max(pc1n);
cols = lines(NoD); %cols = jet(n);

sc_t = [log2(squeeze(mean(cs(18,:,:),2))) pc_p.score(:,2)]; %scatter temp
sc = NaN([size(sc_t) NoD]);

if plt_seg
    figure('Position',[plot_where 800 800]), hold on
    set(gca, 'FontSize', 16)
    set(gcf,'defaultLineLineWidth',2)
    
    xlabel('log2(S/I)')
    ylabel('PC2')
    zlabel('(PC1-min(PC1)).^(1/4)','Interpreter','none')
    axis tight
    cla
    
    for i=1:NoD
        block(i,[1 2]) = [(i-1)*(m1/NoD), i*(m1/NoD)]; %compute lower and upper bounds for block
        sc(and(pc1n>=(block(i,1)),pc1n<block(i,2)),:,i) = sc_t(and(pc1n>=(block(i,1)),pc1n<block(i,2)),:); %sorts values into order based on block membership
        fit_t(i,:) = polyfit(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i),1); %fitty, lol. Seriously though, temp.
        
        x_seg = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
        y_seg = (fit_t(i,1) * x_seg) + fit_t(i,2);
        
        if any(fit_t(i,:))
            scatter3(sc(:,1,i),sc(:,2,i),pc1n,[],cols(i,:),'.');
            plot3(x_seg,y_seg,repmat(mean(block(i,:)),size(x_seg,2),1),'Color',cols(i,:));
        end
        drawnow
        %pause(0.2)
    end
end

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%

plt_ass = 1;

if plt_ass
    figure('Position',[plot_where 800 800]), hold on
    set(gca, 'FontSize', 16)
    set(gcf,'defaultLineLineWidth',2)
    
    x_ass = mean(block(and(fit_t(:,1),fit_t(:,2)),:),2);
    y_1 = fit_t(and(fit_t(:,1),fit_t(:,2)),1);
    y_2 = fit_t(and(fit_t(:,1),fit_t(:,2)),2);
    plot(x_ass,y_1,'r')
    plot(x_ass,y_2,'b')
    
    p_1 = polyfit(log(x_ass),log(y_1),2);
    p_2 = real(polyfit(log(x_ass),log(y_2),2)); %goes imaginary due to negative numbers (I think?)
    y_1p = polyval(p_1, log(x_ass));
    y_2p = polyval(p_2, log(x_ass));
    
    plot(x_ass,exp(y_1p),'r:')
    plot(x_ass,exp(y_2p),'b:')
    
    
    legend('m','c',num2str(p_1),num2str(p_2),'Location','Northwest')
end

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure