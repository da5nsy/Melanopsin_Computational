clc, clear, close all

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
T_SRF = sur_vrhel';
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

%% Modify colour signals
mod_I = 0;

if mod_I
    figure, hist(LMSRI(5,:),50)
    %LMSRI(5,LMSRI(1,:)<0.1) = 0.0000000000000001;
    LMSRI(5,:) = log(LMSRI(5,:))-min(log(LMSRI(5,:)));
    figure, hist(LMSRI(5,:),50)
end

%% PCA of illums

[pc.COEFF, pc.SCORE, pc.LATENT, pc.TSQUARED, pc.EXPLAINED] = pca(T_SPD','VariableWeights',pc.variableweights);

%% Calculate lines of fit between S and I values

clear fit
for i=1:size(LMSRI,3)
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
end %need to visualise this to make sure it's working as planned !!!!!!!!!!

f_m = permute(repmat(fit(:,1),1,1,170),[2,3,1]);
f_c = permute(repmat(fit(:,2),1,1,170),[2,3,1]);

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
    );

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

ds = 1; %downsample

% S/I against PC 2
figure, hold on
scatter3(squeeze(mean(cs(18,:,:),2)),pc.SCORE(:,2),squeeze(mean(LMSRI(5,:,:))),'r.')
set(gca,'Color',repmat(c_norm(18,2),3,1))
axis tight
xlabel(plt_lbls{18})
ylabel('PC 2 weight')
zlabel('I')

figure, hold on
scatter3(squeeze(mean(cs(18,:,:),2)),pc.SCORE(:,2),(pc.SCORE(:,1)-min(pc.SCORE(:,1))).^(1/4),'b.')
set(gca,'Color',repmat(c_norm(18,2),3,1))
axis tight
xlabel(plt_lbls{18})
ylabel('PC 2 weight')
zlabel('PC 1 weight')
axis equal

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
for j=50:2:80
    
    msi = squeeze(mean(cs(18,:,:),2))-0.5867; %mean S/I
    s = squeeze(mean(cs(2,:,:),2));
    pc2 = pc.SCORE(:,2);
    Z = (pc.SCORE(:,1)-min(pc.SCORE(:,1))).^(1/4);
    
    r=zeros(2,2600);
    
    for i=1:2600
        theta = (8 - Z(i))*j;
        R = [cosd(theta) -sind(theta); ...
            sind(theta) cosd(theta)];
        r(:,i) = R*[msi(i) pc2(i)]';
    end
    
    figure,
    subplot(1,2,1)
    scatter3(msi,pc2,Z,'k.')
    view(2)
    axis equal
    
    subplot(1,2,2)
    scatter3(r(1,:),r(2,:),Z,'k.')
    view(2)
    axis equal
    
    title(j)
    
    pause(1)
end
