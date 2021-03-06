function [corr_return, p_1, p_2] = melcomp_3(fv_ind,sv_ind)

% fv_ind = 18; % first value index
% sv_ind = 1; %second value index

clear, clc, close all

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_3';

set(groot,'defaultfigureposition',[100 100 500 400]); 
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesFontName', 'Courier');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vectors
set(groot,'defaultfigurecolor','white')

ff = '-dpng'; %file format
p  = 0; %print? (aka save?), set to 1 to begin saving

plot_where = [800,50];
plot_size  = [800,375];

%% Start the figure

figure
hold on
set(gca, 'FontSize', 16)

%% Load data

[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader('SRF','Vrhel_nat_extended','SSF','SP','Lum','SP');

% Compute colorimetry
[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

% Log conversion straightens out some of the relationships for easier
% handling later
lsri = log(lsri);

%% SPDs (Spectral Power Distributions)

% try load melcomp_3_addI.mat %bit slow so I saved it out
% catch
%     % Additional SPD info (source: personal correspondance with J. Hern�ndez-Andr�s)
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
    axis tight
    xticks(400:100:700);
    ylim([0,2.1]),yticks([0,1,2]);
    xlabel('Wavelength (nm)'); ylabel('Power (W m ^-^2 nm^-^1)');
    
    plot(SToWls(S_sh),T_SPD(:,103)) %single spd
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),T_SPD) %all spd
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%% SRFs (Spectral Reflectance Functions)

plt_refs = 1;
if plt_refs
    cla
    xticks(400:100:700);
    ylim([0 1]); yticks([0,1]);
    xlabel('Wavelength (nm)'); ylabel('Reflectance');
    
    plot(SToWls(S_sh),T_SRF(:,72)) %pear
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),T_SRF) %all
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

% T_SRF = mean(T_SRF,2);
%This is to check the effect of averaging the
%reflectances at this stage rather than later in the process.
%Spolier: nothing changes! (except a couple of errors being thrown.

%% SSF (Spectral Sensitivity Functions)

plt_spec_sens = 1;
if plt_spec_sens
    figure(1)
    cla
    xlim([390 730]); xticks(400:100:700);
    ylim([0 1]); yticks([0,1]);
    ylabel('Sensitivity')
    plot(SToWls(S_sh),T_SSF)
    %plot(SToWls(S_rods),T_rods)
    %plot(SToWls(S_mel),T_mel)
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%% compute peaks of spectral sensitivities
[~,ms_ind] = max(T_SSF);
S_full = SToWls(S_sh);
ms_val = S_full(ms_ind);

%% PCA

% PCA of spectral Power distributions
% compute pca variable weight
vw = 0;
if vw == 0 %no variable weights
    pc_p.variableweights = ones(81,1);
elseif vw == 1 %S + L and max in between
    pc_p.variableweights = T_SSF(:,1)+T_SSF(:,3);
    [~,t_p_locs] = findpeaks(pc_p.variableweights);
    pc_p.variableweights(t_p_locs(1):t_p_locs(2)) = max(pc_p.variableweights);
    pc_p.variableweights(pc_p.variableweights>1) = 1;
elseif vw == 2 %S+L+I
    pc_p.variableweights = T_SSF(:,1)+T_SSF(:,3)+T_mel(:,1);
end

plt_vw = 1;
if plt_vw
    cla
    plot(SToWls(S_sh),pc_p.variableweights)
    xlim([390 730]);
    xticks(400:100:700);
    ylim([0 1]);
    yticks(ylim);
    ylabel('Weight')    
    if p, print([base,'\vw'],ff); p=p+1; end %save figure
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
        'k','LineStyle','none','FaceAlpha','0.03','HandleVisibility','off');
    
    %normalised PCs:
    %plot(SToWls(S_sh),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))
    
    plot(SToWls(S_sh),pc_p.coeff(:,1),'DisplayName','PC1')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),pc_p.coeff(:,2),'DisplayName','PC2')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),pc_p.coeff(:,3),'DisplayName','PC3')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    hold on
    ylim([0 2]); yticks([0,1,2]);
    
    %reconstructed data
    recon_p(:,1) = pc_p.score(103,:)   * pc_p.coeff'        + pc_p.mu; % Reconstructed with all PCs
    recon_p(:,2) = pc_p.score(103,1:3) * pc_p.coeff(:,1:3)' + pc_p.mu; % Reconstructed with 3 PCs
    
    %plot(SToWls(S_sh), T_SPD(:,103), 'b') % Original data
    %plot(SToWls(S_sh), recon_p(:,1), 'g') % All PC
    plot(SToWls(S_sh), recon_p(:,2), 'r') % 3 PC
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
    
    %plot(SToWls(S_sh),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))
    
    plot(SToWls(S_sh),pc_r.coeff(:,1))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),pc_r.coeff(:,2))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot(SToWls(S_sh),pc_r.coeff(:,3))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    hold on
    ylim([0 1]); yticks([0,1]);
    
    recon_r(:,1) = pc_r.score(72,:)   * pc_r.coeff'        + pc_r.mu; % Reconstructed with all PCs
    recon_r(:,2) = pc_r.score(72,1:3) * pc_r.coeff(:,1:3)' + pc_r.mu; % Reconstructed with 3 PCs
    
    %plot(SToWls(S_sh), T_SRF(:,72), 'b') % Original data
    %plot(SToWls(S_sh), recon_r(:,1), 'g') % All PC
    plot(SToWls(S_sh), recon_r(:,2), 'r') % 3 PC
    % legend({'Original Data','Reconstructed with 1:3'})
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    % figure, plot((pc_r.score*pc_r.coeff' + repmat(pc_r.mu,length(pc_r.score),1))') %all the data reconstructed
    
    % disp(pc_r.score(72,1:3)) %values for demo ref
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
    
    plot(SToWls(S_sh),pc_p.coeff(:,2),'Color',cols(2,:))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot([ms_val(3),ms_val(3)],[min(ylim),max(ylim)],'Color','k')
    plot([ms_val(5),ms_val(5)],[min(ylim),max(ylim)],'Color','k')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    plot(SToWls(S_sh),pc_p.coeff(:,2)*-1,'Color',cols(2,:)*((-1+2)/4))
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    plot([ms_val(3),ms_val(3)],[min(ylim),max(ylim)],'Color','k')
    plot([ms_val(5),ms_val(5)],[min(ylim),max(ylim)],'Color','k')
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
    
    cla
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    for i = -1.5:0.5:1.5
        plot(SToWls(S_sh),pc_p.coeff(:,2)*i,'Color',cols(2,:)*((i+3)/4))
    end
    if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end


%% Calculate lines of fit between S and I values

for i=1:size(LMSRI,3)
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
end

f_m = permute(repmat(fit(:,1),1,1,size(LMSRI,2)),[2,3,1]); % Mx + c
f_c = permute(repmat(fit(:,2),1,1,size(LMSRI,2)),[2,3,1]); % mx + C

%% - % Checking

% figure, hold on
% scatter(log10(LMSRI(3,:,1)),log10(LMSRI(5,:,1)))
% x = linspace(-2,0);
% y = fit(1,1)*x;
% plot(x,y)
%
%
% figure, hold on
% scatter(log10(LMSRI(3,:,2)),log10(LMSRI(5,:,2)))
% x = linspace(-2,0);
% y = repmat(fit(2,2), length(x),1);
% plot(x,y)
%
% % - %
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
%      x_t2 = linspace(min(log10(LMSRI(3,:,i))),max(log10(LMSRI(3,:,i))));
%      y_t2 = x_t*f_m(1,1,i) + f_c(1,1,i);
%
%     %plot3(x_t,y_t,ones(100,1)*i,'r')
%     plot3(x_t2,y_t2,ones(100,1)*i,'g:')%
% end
%
% %view(-43,0)
% %
% % figure,
% % scatter3(f_m(1,1,:),f_c(1,1,:),pc_p.score(:,2))
% % xlabel('f_m'); ylabel('f_c');

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
plt_lbls{11} = 'SP\_Lum'; %'(0.6373*L)+(0.3924*M)';
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

exploratory_figs = 0;

if exploratory_figs %open question of whether z should be pc1 or L
    
    % mean of data (as above):
    figure,
    scatter3(squeeze(mean(cs(18,:,:),2)),...
        pc_p.score(:,2),...
        pc_p.score(:,1),...
        'k.')
    
    % all data (for S/I, aka 18)
    figure, hold on
    scatter3(cs(18,:),...
        repelem(pc_p.score(:,2),size(cs,2)),...
        repelem(pc_p.score(:,1),size(cs,2)),...
        'k.')
    
    % all data for all signals
    for fv_ind=1:size(cs,1)
        figure,
        scatter3(cs(fv_ind,:),...
            repelem(pc_p.score(:,2),size(cs,2)),...
            repelem(pc_p.score(:,1),size(cs,2)),...
            'k.')
        
        xlabel(plt_lbls(fv_ind))
        ylabel('PC2')
        zlabel('PC1')
    end
    
    % all data for S/I removing low values
    
    removeLowVals = 1;
    
    if removeLowVals
        x = cs(18,cs(1,:) > 0.5); %the 18th signal (S/I) where the first signal (L) is greater than 0.5
        y1 = repelem(pc_p.score(:,2),size(cs,2))';
        y = y1(cs(1,:) > 0.5);
        z = cs(1,cs(1,:) > 0.5);
    else
        x = cs(18,:); %the 18th signal (S/I) where the first signal (L) is greater than 0.5
        y = repelem(pc_p.score(:,2),size(cs,2))';
        
        z = cs(1,:);
    end
    
    figure,
    scatter3(x,y,z,'k.')
    %scatter3(x,y,z,10,'MarkerEdgeColor',[0.1,0.1,0.1],'MarkerEdgeAlpha',0.5)
    
    xlabel(plt_lbls(18))
    ylabel('PC2')
    zlabel(plt_lbls(1))
    
end

crl = zeros(size(cs,1),nPC);
for j=1:nPC
    for i=1:size(cs,1)
        crl(i,j) = corr(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j));
    end
end
crl=abs(crl);

crl_norm = crl;
for i=1:j %leftover, be careful
    crl_norm(:,i) = crl_norm(:,i) - min(crl_norm(:,i));
    crl_norm(:,i) = crl_norm(:,i)/max(crl_norm(:,i));
end

plt_viz = 1;
if plt_viz
    figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
    hold on
    
    subplot(1,2,1)
    imagesc(crl)
    colorbar
    title({'correlation between signal', 'and PC weight'})
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:size(cs,1));
    set(gca, 'YTickLabel', plt_lbls);
    colormap('gray');
    %set(gca, 'FontSize', 16)
    
    subplot(1,2,2)
    imagesc(crl_norm)
    colorbar
    title('normalised')
    xlabel('PC')
    
%     set(gca, 'XTick', 1:nPC);
     set(gca, 'YTick', []);
%     set(gca, 'YTickLabel', plt_lbls);
%     colormap('gray');
%     set(gca, 'FontSize', 16)
    
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
            set(gca,'Color',repmat(crl_norm(i,j),3,1))
            axis tight
        end
    end
end


plt_sca2 = 1;
if plt_sca2
    counter = 1;
    figure('Position',[plot_where 250 900])
    for i=[6:9, 13:18]
        for j=2
            subplot(length([6:9, 13:18]),1,counter), hold on
            counter = counter + 1;
            scatter(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j),'r.')
            %ylabel(plt_lbls{i},'rotation',0)
            legend(plt_lbls{i},'AutoUpdate', 'Off')
            set(gca,'XTick',[])
            set(gca,'XTickLabel',[])
            set(gca,'YTick',[])
            set(gca,'YTickLabel',[])
            
            set(gca,'Color',repmat(crl_norm(i,j),3,1))
            axis tight
            
            pft = polyfit(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j),1);
            x = linspace(min(squeeze(mean(cs(i,:,:),2))),max(squeeze(mean(cs(i,:,:),2))));
            y = pft(1)*x +pft(2);
            plot(x,y,'k')
        end
    end
end

%% Scatter plots between...

%L against PC1
figure(1), hold on
cla
scatter(squeeze(mean(cs(1,:,:),2)),pc_p.score(:,1),'k.')
xlabel(plt_lbls{1})
ylabel('PC1')
xlim('auto')
xticks('auto')
ylim('auto')
yticks('auto')

f_L1 = polyfit(squeeze(mean(cs(1,:,:),2)),pc_p.score(:,1),1);
x = linspace(min(squeeze(mean(cs(1,:,:),2))),max(squeeze(mean(cs(1,:,:),2))));
y = x*f_L1(1)+f_L1(2);
plot(x,y,'b')

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

% S/I against PC2
cla
scatter(squeeze(mean(cs(18,:,:),2)),pc_p.score(:,2),'k.')
xlabel(plt_lbls{18})
ylabel('PC2')

f_SI2 = polyfit(squeeze(mean(cs(18,:,:),2)),pc_p.score(:,2),1);
x = linspace(min(squeeze(mean(cs(18,:,:),2))),max(squeeze(mean(cs(18,:,:),2))));
y = x*f_SI2(1)+f_SI2(2);
plot(x,y,'b')

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% 3D plot - basic (same as previous but without line and with a third dimension)

if and(exist('fv_ind','var'),exist('sv_ind','var'))
else
    fv_ind = 18; % first value index
    sv_ind = 1; %second value index
    disp(['using default values for fv_ind and sv_ind: ',num2str(fv_ind),', ',num2str(sv_ind)])
end

fv = squeeze(mean(cs(fv_ind,:,:),2));
sv = squeeze(mean(cs(sv_ind,:,:),2));

cla
scatter3(fv,pc_p.score(:,2),sv,'k.')

xlabel(plt_lbls{fv_ind})
ylabel('PC2')
zlabel(plt_lbls{sv_ind})

% needs to be 3D, make gif? !!!!!!!!!!!!
% if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Plot showing segmentation by sv value

plt_seg = 1;

NoD = 40; %nnumber of divisions
mI = max(sv);
cols = lines(NoD); %cols = jet(n);

sc_t = [fv pc_p.score(:,2)]; %scatter temp
sc = NaN([size(sc_t) NoD]);

if plt_seg
    figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
    hold on
    set(gca, 'FontSize', 16)
    
    xlabel(plt_lbls{fv_ind})
    ylabel('PC2')
    zlabel(plt_lbls{sv_ind})
    
    axis tight
    cla
    
    for i=1:NoD
        block(i,[1 2]) = [(i-1)*(mI/NoD), i*(mI/NoD)]; %compute lower and upper bounds for block
        sc(and(sv>=(block(i,1)),sv<block(i,2)),:,i) = sc_t(and(sv>=(block(i,1)),sv<block(i,2)),:); %sorts values into order based on block membership
        fit_t(i,:) = polyfit(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i),1); %fit temp
        
        x_seg = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
        y_seg = (fit_t(i,1) * x_seg) + fit_t(i,2);
        
        if any(fit_t(i,:))
            scatter3(sc(:,1,i),sc(:,2,i),sv,[],cols(i,:),'.');
            plot3(x_seg,y_seg,repmat(mean(block(i,:)),size(x_seg,2),1),'Color',cols(i,:));
        end
        drawnow
        %pause(0.2)
    end
end

% would be better as 3D gif !!!!!!!!!!!
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Assess trends in fitted lines through data
plt_ass = 1;

if plt_ass
    figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
    hold on
    set(gca, 'FontSize', 16)
    
    fit_t_idx = and(fit_t(:,1),fit_t(:,2));
    x_ass = mean(block(fit_t_idx,:),2);
    y_1 = fit_t(fit_t_idx,1);
    y_2 = fit_t(fit_t_idx,2);
    scatter(x_ass,y_1,'ro')
    scatter(x_ass,y_2,'bo')
    
    xlabel(plt_lbls{sv_ind})
    ylabel('Value in line of best fit')
    
    p_1 = polyfit(x_ass,y_1,1);
    p_2 = polyfit(x_ass,y_2,1);
    y_1p = polyval(p_1, x_ass);
    y_2p = polyval(p_2, x_ass);
    
    plot(x_ass,y_1p,'r:')
    plot(x_ass,y_2p,'b:')
    
    
    legend('m','c',num2str(p_1),num2str(p_2),'Location','Northwest')
end

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%% Visualise model fit to data

figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
hold on
set(gca, 'FontSize', 16)

scatter3(fv,pc_p.score(:,2),sv,'k.')
xlabel(plt_lbls{fv_ind})
% ylabel({'PC2',...
%     ['or (',num2str(p_1(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_1(2),3),') * ',plt_lbls{fv_ind},...
%     ' + ','(',num2str(p_2(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_2(2),3),')']})
zlabel(plt_lbls{sv_ind})

for i=linspace(min(sv),max(sv))
    x = linspace(min(fv),max(fv));
    m = (p_1(1) * i) + p_1(2); %fits better without second terms but there's no logical argument for why to exclude them
    c = (p_2(1) * i) + p_2(2);
    y = m*x + c;
    plot3(x,y,ones(100,1)*i,'r')
end

% needs to be 3D, make gif? !!!!!!!!!!!!
% if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

legend({'Data: PC2',['Model: (',num2str(p_1(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_1(2),3),') * ',plt_lbls{fv_ind},...
    ' + (',num2str(p_2(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_2(2),3),')']})

% %hard code print with additional prefix
%print([base,'\f6\f6_',num2str(fv_ind),'_',num2str(sv_ind)],ff);

%% Consider 2D correlation of model to PC2

figure
hold on
set(gca, 'FontSize', 16)

m = p_1(1) * sv + p_1(2);
c = p_2(1) * sv + p_2(2);

estimatedPC2 =  m .* (fv) + c;

%scatter(estimatedPC2,pc_p.score(:,2),'k.')

%scatter3(estimatedPC2,pc_p.score(:,2),sv,'k.')
%plot(estimatedPC2,polyval(polyfit(estimatedPC2,pc_p.score(:,2),1),estimatedPC2),'k')

scatter3(estimatedPC2(sv>0.5),pc_p.score((sv>0.5),2),sv(sv>0.5),'k.')
y = polyval(polyfit(estimatedPC2(sv>0.5),pc_p.score((sv>0.5),2),1),estimatedPC2(sv>0.5));
plot(estimatedPC2(sv>0.5),y,'b')

disp(polyfit(estimatedPC2(sv>0.5),pc_p.score((sv>0.5),2),1))

%xlabel(['Estimated PC2 based on ' plt_lbls{fv_ind}])
xlabel({['(',num2str(p_1(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_1(2),3),') * ',plt_lbls{fv_ind}],...
    [' + (',num2str(p_2(1),3),' * ',plt_lbls{sv_ind},' + ',num2str(p_2(2),3),')']})
ylabel('PC2')
zlabel(plt_lbls{sv_ind})
cleanTicks

% %hard code print with additional prefix
%print([base,'\f7\f7_',num2str(fv_ind),'_',num2str(sv_ind)],ff);

if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

%%
corr_return = corr(estimatedPC2(sv>0.5),pc_p.score((sv>0.5),2));

return

% The following currently crashes MATLAB

% %% Scatter FULL data
% 
% figure, hold on
% for ref = 1:size(T_SRF,2)
%     scatter3(cs(18,ref,:),pc_p.score(:,2),pc_p.score(:,1),'.')
% end
% 
% figure, hold on
% for ref = 1:size(T_SRF,2)
%     scatter3(cs(14,ref,:),pc_p.score(:,2),pc_p.score(:,1),'.')
% end
% 
%% New version of 'Calc correlation between' taking into account amended signals

clear, clc, close all
lore = load('melcomp_3_correlation_results.mat'); %loaded results

plt_viz = 1;
if plt_viz
    figure
    hold on
    
    ncs = size(lore.corr_return,1);%number of colour signals
    nPC = 3; %number of principal components
    crl2 = zeros(ncs,nPC);
    crl2(:,1:3) = lore.corr_return(:,1:3);
    
    subplot(1,2,1)
    imagesc(crl2)
    colorbar
    title('correlation between signal and PC weight')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:ncs);
    set(gca, 'YTickLabel', lore.plt_lbls);
    colormap('gray');
    
    crl_norm2 = crl2;
    for i=1:3 %leftover, be careful
        crl_norm2(:,i) = crl_norm2(:,i) - min(crl_norm2(:,i));
        crl_norm2(:,i) = crl_norm2(:,i) / max(crl_norm2(:,i));
    end
    
    subplot(1,2,2)
    imagesc(crl_norm2)
    colorbar
    title('normalised')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:ncs);
    set(gca, 'YTickLabel', lore.plt_lbls);
    colormap('gray');
    
    %if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end

%%
plt_sca2 = 1;
if plt_sca2
    counter = 1;
    figure('Position',[plot_where 250 900])
    for i=[6:9, 13:18]
        for j=2
            subplot(length([6:9, 13:18]),1,counter), hold on
            counter = counter + 1;
            %scatter(squeeze(cs(i,:,:)),repmat(pc_p.score(:,j),1,83,1)','r.')
            scatter(reshape(squeeze(cs(i,:,:)),1,[]),reshape(repmat(pc_p.score(:,j),1,83,1),1,[])','r.')
            %ylabel(plt_lbls{i},'rotation',0)
            legend(plt_lbls{i},'AutoUpdate', 'Off')
            set(gca,'XTick',[])
            set(gca,'XTickLabel',[])
            set(gca,'YTick',[])
            set(gca,'YTickLabel',[])
            
            set(gca,'Color',repmat(crl_norm(i,j),3,1))
            axis tight
            
            pft = polyfit(squeeze(mean(cs(i,:,:),2)),pc_p.score(:,j),1);
            x = linspace(min(squeeze(mean(cs(i,:,:),2))),max(squeeze(mean(cs(i,:,:),2))));
            y = pft(1)*x +pft(2);
            plot(x,y,'k')
        end
    end
end


end
