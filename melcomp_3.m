function [cs,pc,plt_lbls] = melcomp_3

try
    nargin;
catch
    clear, clc, close all
end

%% SPDs (Spectral Power Distributions)

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)
T_SPD=final; clear final
S_SPD=[300,5,161];

try load SPD.mat addI %bit slow so I saved it out
catch
    % Additional SPD info (source: personal correspondance with J. Hernández-Andrés)
    [addI.NUM,addI.TXT,addI.RAW] = xlsread('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\add_info.xlsx');
    for i=1:length(T_SPD) %a lot of stupid code just to get the date and time out
        addI.t(i) = datetime(...
            [char(addI.RAW(2,2)),' ',char(days(cell2mat(addI.RAW(2,3))),'hh:mm')],...
            'InputFormat','dd/MM/uuuu HH:mm');
    end
    addI.el = addI.NUM(:,4); %elevation
    addI.az = addI.NUM(:,5); %azimuth
end

% T_SPD = T_SPD(:,addI.el>30);
% addI.el = addI.el(addI.el>30);


%% SRFs (Spectral Reflectance Functions)

load sur_vrhel
T_SRF = sur_vrhel(:,[1:44,65,69,81:154])';
S_SRF = S_vrhel;
clear sur_vrhel S_vrhel

%% SSF (Spectral Sensitivity Functions)

% change for ss at some point in future !!!!!!!!!!1
load T_cones_sp T_cones_sp S_cones_sp
T_SSF = T_cones_sp;
S_SSF = S_cones_sp;
clear T_cones_sp S_cones_sp

load T_rods T_rods S_rods
load T_melanopsin T_melanopsin S_melanopsin
T_mel = SplineCmf(S_melanopsin,T_melanopsin, S_melanopsin - [10, 0, 0],1); %Increasing the range of this function in case it ends up limiting the range of S_sh, and shorten variable names
S_mel = S_melanopsin - [10, 0, 0]; clear S_melanopsin T_melanopsin

%% Interpolate functions to match range and interval

%reduce all data down to the common range/interval
S_sh = [max([S_SPD(1),S_SRF(1),S_SSF(1)]),max([S_SPD(2),S_SRF(2),S_SSF(2)]),min([S_SPD(3),S_SRF(3),S_SSF(3)])]; %S_shared: work out what the lowest common denominator for the range/interval of the data is
%S_sh = [400,5,61]; % brute force way to do something like the variableweights thing

T_SPD = SplineSpd(S_SPD,T_SPD,S_sh,1); % extend == 1: Cubic spline, extends with last value in that direction
T_SRF = SplineSrf(S_SRF,T_SRF',S_sh,1);
T_SSF  = SplineCmf(S_SSF,T_SSF,S_sh,1)';
T_rods = SplineCmf(S_rods,T_rods,S_sh)'; %extend? !!!!!!!!!!
T_mel  = SplineCmf(S_mel,T_mel,S_sh)'; %extend? !!!!!!!!!!
[S_SPD, S_SRF, S_SSF, S_rods, S_mel] = deal(S_sh);

% combine sensitivity vectors
T_LMSRI=[T_SSF,T_rods,T_mel];
S_LMSRI=S_sh;

%% compute pca variable weight
vw = 1;

if vw == 0 %no variable weights
    pc.variableweights = ones(81,1);
elseif vw == 1 %S + L and max in between
    pc.variableweights = T_SSF(:,1)+T_SSF(:,3); 
    [~,t_p_locs] = findpeaks(pc.variableweights);
    pc.variableweights(t_p_locs(1):t_p_locs(2)) = max(pc.variableweights);
    %plot(pc.variableweights)
elseif vw == 2 %S+L+I
    pc.variableweights = T_SSF(:,1)+T_SSF(:,3)+T_mel(:,1); 
end

%% Calculate colour signals

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_LMSRI'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i)); %change this to be more direct !!!!!!!!!!!!!!
    lsri(3,:,i)   = LMSRI(4,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); %change this to be more direct !!!!!!!!!!!!!!
    lsri(4,:,i)   = LMSRI(5,:,i)./(0.6373*LMSRI(1,:,i)+0.3924*LMSRI(2,:,i)); %change this to be more direct !!!!!!!!!!!!!!
end

%% PCA of illums

[pc.COEFF, pc.SCORE, pc.LATENT, pc.TSQUARED, pc.EXPLAINED] = pca(T_SPD','VariableWeights',pc.variableweights);

%% Calculate lines of fit between S and I values

clear fit
for i=1:size(LMSRI,3)
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
end %need to visualise this to make sure it's working as planned !!!!!!!!!!

f_m = permute(repmat(fit(:,1),1,1,size(LMSRI,2)),[2,3,1]);
f_c = permute(repmat(fit(:,2),1,1,size(LMSRI,2)),[2,3,1]);

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
        c(i,j) = corr(squeeze(mean(cs(i,:,:),2)),pc.SCORE(:,j));
    end
end
c=abs(c);

c_norm = c;
for i=1:j %leftover, be careful
    c_norm(:,i) = c_norm(:,i) - min(c_norm(:,i));
    c_norm(:,i) = c_norm(:,i)/max(c_norm(:,i));
end

%%
plt_viz = 1;

if plt_viz
    figure, hold on
    subplot(1,2,1)
    imagesc(c)
    colorbar
    title('correlation between signal and PC weight')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:size(cs,1));
    set(gca, 'YTickLabel', plt_lbls);
    colormap('gray');
    
    subplot(1,2,2)
    imagesc(c_norm)
    colorbar
    title('as above, normalised')
    xlabel('PC')
    
    set(gca, 'XTick', 1:nPC);
    set(gca, 'YTick', 1:size(cs,1));
    set(gca, 'YTickLabel', plt_lbls);
    colormap('gray');
end

%% plot scatters with backrounds the colours of the strength of correlation

plt_sca = 0;

if plt_sca
    figure,
    for i=1:size(cs,1)
        for j=1:nPC
            subplot(size(cs,1),nPC,j+((i-1)*nPC))
            scatter(squeeze(mean(cs(i,:,:),2)),pc.SCORE(:,j),'r.')
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

%% plot individual scatters of interest

% ds = 1; %downsample

% % S/I against PC 2
% scatter(squeeze(mean(cs(1,:,:),2)),pc.SCORE(:,1),'k.')
% xlabel(plt_lbls{1})
% ylabel('PC 1')
% 
% scatter3(squeeze(mean(cs(18,:,:),2)),pc.SCORE(:,2),pc.SCORE(:,1),'k.')
% %set(gca,'Color',repmat(c_norm(18,2),3,1))
% xlabel(plt_lbls{18})
% ylabel('PC 2')
% zlabel('PC 1')

% figure, hold on
% scatter3(squeeze(mean(cs(18,:,:),2)),pc.SCORE(:,2),(pc.SCORE(:,1)-min(pc.SCORE(:,1))).^(1/4),'b.')
% set(gca,'Color',repmat(c_norm(18,2),3,1))
% axis tight
% xlabel(plt_lbls{18})
% ylabel('PC 2 weight')
% zlabel('PC 1 weight')
% axis equal

% yfit = polyval(mean(fit),log10(LMSRI(3,:,i)));
% plot(log10(LMSRI(3,:,i)),yfit,'Color',cols(i,:))

% % l against PC3
% figure,
% scatter3(squeeze(mean(cs(6,:,:),2)),pc.SCORE(:,3),addI.el,'r.')
% set(gca,'Color',repmat(c_norm(6,3),3,1))
% axis tight
% xlabel(plt_lbls{6})
% ylabel('PC 3 weight')

% % c against PC2 (logs)
% figure,
% scatter3(squeeze(mean(cs(21,:,:),2)),log10(pc.SCORE(:,2)),log10(addI.el),'r.')
% set(gca,'Color',repmat(c_norm(21,2),3,1))
% axis tight
% xlabel(plt_lbls{21})
% ylabel('PC 2 weight')

% figure,
% %plot(pc.SCORE(1:ds:end,1))
% scatter(addI.el(1:ds:end,1),pc.SCORE(1:ds:end,2))
% xlabel('Elevation')
% ylabel('PC2')

%%

% figure,
% plot(SToWls(S_SPD),pc.COEFF(:,1:3)./max(pc.COEFF(:,1:3)))
% legend({'PC1','PC2','PC3'},'Location','Best')
% axis tight

%%
plt_soil = 0;

SoI = squeeze(mean(cs(18,:,:),2)); %signal of interest, OR S-over-I (take your pick!)
SoIL = log2(SoI); %Log SoI
pc1 = pc.SCORE(:,1);
pc2 = pc.SCORE(:,2);
pc1n = (pc1-min(pc1)).^(1/4); %pc1 normalised

% figure, hold on
% scatter3(SoI,pc2,pc1,'r.')
% axis tight
% xlabel('SoI')
% ylabel('PC2')
% zlabel('PC1')
%
% figure, hold on
% scatter3(SoI,pc2,pc1n,'b.')
% axis tight
% xlabel('SoI')
% ylabel('PC2')
% zlabel('PC1n')
% axis equal

if plt_soil
    figure, hold on
    scatter3(SoIL,pc2,pc1n,'b.')
    axis tight
    xlabel('SoIL')
    ylabel('PC2')
    zlabel('PC1n')
    axis equal
end
%% Segment

plt_seg = 0;

NoD = 30; %nnumber of divisions
m1 = max(pc1n);
cols = flag(NoD); %cols = jet(n);

sc_t = [SoIL pc2]; %scatter temp
sc = NaN([size(sc_t) NoD]);

if plt_seg
    figure, hold on
    axis equal
    
    for i=1:NoD
        sc(and(pc1n<(i*(m1/NoD)),pc1n>((i-1)*(m1/NoD))),:,i) = sc_t(and(pc1n<(i*(m1/NoD)),pc1n>((i-1)*(m1/NoD))),:);
        fit_t(i,:) = polyfit(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i),1); %fitty, lol. Seriously though, temp.
        
        x_seg = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
        y = (fit_t(i,1) * x_seg) + fit_t(i,2);
        
        if any(fit_t(i,:))
            scatter(sc(:,1,i),sc(:,2,i),[],cols(i,:));
            plot(x_seg,y,'Color',cols(i,:));
        end
        drawnow
        %pause(0.2)
    end
end

%% Assess curve of fit values
plt_ass = 0;

if plt_ass
    figure, hold on
    
    x_ass = 4:30;
    y_1 = fit_t(4:30,1);
    y_2 = fit_t(4:30,2);
    plot(x_ass,y_1,x_ass,y_2)
    
    p_1 = polyfit(log(x_ass),log(y_1'),4);
    p_2 = real(polyfit(log(x_ass),log(y_2'),4)); %goes imaginary due to negative numbers (I think?)
    y_1p = polyval(p_1, log(x_ass));
    y_2p = polyval(p_2, log(x_ass));
    
    plot(x_ass,exp(y_1p),x_ass,exp(y_2p))
    
    legend
end

%plot(x,exp(p(2))*x.^p(1));

% f = fit(x,y,'power2') %need curve fitting toolbox for this


end