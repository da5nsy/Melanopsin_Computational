function melcomp_why(plt)
try
    nargin;
catch
    clear, clc, close all
end

load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141]; %natural ones
T_refs_all = sur_vrhel';
T_refs_nat = sur_vrhel(:,refs)';
S_refs = S_vrhel; clear S_vrhel

load('C:\Users\cege-user\Dropbox\UCL\Data\Reference Data\Granada Data\Granada_daylight_2600_161.mat');
% http://colorimaginglab.ugr.es/pages/Data#__doku_granada_daylight_spectral_database
% From: J. Hernández-Andrés, J. Romero& R.L. Lee, Jr., "Colorimetric and
%       spectroradiometric characteristics of narrow-field-of-view
%       clear skylight in Granada, Spain" (2001)

% T_SPD=final; clear final
% S_SPD=[300,5,161];

%match to cones instead
T_SPD = final(17:97,:); clear final;
S_SPD = [380,5,81];

% figure, hold on;
% plot(SToWls(S_SPD),T_SPD)
% %drawnow, pause(0.3)
% xlabel('Wavelength (nm)')
% ylabel('Spectral Power Distribution (W m ^-^2 nm^-^1)')
% xlim([S_SPD(1),max(SToWls(S_SPD))]), ylim([0 max(T_SPD(:))]);

%normalise T_SPD for luminance
load T_cones_sp

SPD_L = (0.6373*T_cones_sp(1,:)*T_SPD)+(0.3924*T_cones_sp(2,:)*T_SPD);
T_SPD_n = T_SPD./SPD_L;

% figure, hold on;
% plot(SToWls(S_SPD),T_SPD_n)


%% Testing STD of refs

% % Natural, wavelength
% figure, hold on
% plot(SToWls(S_refs),T_refs_nat)
% plot(SToWls(S_refs),std(T_refs_nat),'ko','DisplayName','STD')
%
% % Natural, frequency
% figure, hold on
% plot(1./SToWls(S_refs),T_refs_nat)
% plot(1./SToWls(S_refs),std(T_refs_nat),'ko','DisplayName','STD')
%
% % % All, wavelength
% % figure, hold on
% % plot(SToWls(S_refs),T_refs_all)
% % plot(SToWls(S_refs),std(T_refs_all),'ko','DisplayName','STD')
% %
% % % All, frequency
% % figure, hold on
% % plot(1./SToWls(S_refs),T_refs_all)
% % plot(1./SToWls(S_refs),std(T_refs_all),'ko','DisplayName','STD')

% %% Testing STD of illums
% % Natural, wavelength
% figure, hold on
% plot(SToWls(S_SPD),T_SPD_n)
% plot(SToWls(S_SPD),std(T_SPD_n'),'ko','DisplayName','STD')
%
% % Natural, frequency
% figure, hold on
% plot(1./SToWls(S_SPD),T_SPD_n)
% plot(1./SToWls(S_SPD),std(T_SPD_n'),'ko','DisplayName','STD')

%% Plot correlation

plt_correl = 0;
try strcmp(plt,'correl');
    if strcmp(plt,'correl')
        plt_correl = 1;
    end
catch
end

load T_cones_ss2.mat T_cones_ss2 S_cones_ss2
load T_melanopsin.mat T_melanopsin S_melanopsin
load T_rods.mat T_rods S_rods
T_rods = SplineCmf(S_rods,T_rods,S_cones_ss2); S_rods = S_cones_ss2;

[~, maxloc]        = max(T_cones_ss2,[],2);
[~, maxloc(end+1)] = max(T_rods,[],2);
[~, maxloc(end+1)] = max(T_melanopsin,[],2);
S_cones_ss2_f = SToWls(S_cones_ss2);
maxl = S_cones_ss2_f(maxloc); %maxl for max lambda

c = corr(T_refs_nat);

if plt_correl
    figure, hold on
    imagesc(c)
    axis image
    colorbar
    
    %fill([0, max(xlim), 0, 0],[0, max(ylim), max(ylim), 0],[1,1,1],'LineStyle','none')
    
    S_refs_f = SToWls(S_refs); %f for full
    set(gca,'XTickLabel',S_refs_f(xticks+1))
    set(gca,'YTickLabel',S_refs_f(xticks+1))
    
    txt='LMSRI';
    for i=0:4
        scatter((maxl-S_cones_ss2(1))/2,(circshift(maxl,i)-S_cones_ss2(1))/2,'k^')
    end
    for i=1:5
        text((maxl(i)-S_cones_ss2(1))/2,(maxl(3)-S_cones_ss2(1))/2-5,txt(i))
        text((maxl(3)-S_cones_ss2(1))/2-10,(maxl(i)-S_cones_ss2(1))/2,txt(i))
    end
    
    xlabel('Wavelength (nm)')
    ylabel('Wavelength (nm)')
end


%% Daylight PCA

plt_PCA = 0;
try strcmp(plt,'PCA');
    if strcmp(plt,'PCA')
        plt_PCA = 1;
    end
catch
end

if plt_PCA
    coeff = pca(T_SPD');
    
    figure, hold on
    linetype={'k-','k--','k:'};
    for i=1:2
        plot(SToWls(S_SPD),coeff(:,i),linetype{i},'LineWidth',1)
    end
    for i=1:5
        plot([maxl(i),maxl(i)],[min(ylim),max(ylim)],'Color',[0.5,0.7,1])
    end
    for i=1:2 %replot just to get on top
        plot(SToWls(S_SPD),coeff(:,i),linetype{i},'LineWidth',1)
    end
    xlim([min(SToWls(S_SPD)),max(SToWls(S_SPD))])
    %plot([min(xlim),max(xlim)],[0,0],'b--') %Plot zero line
    
    xlabel('Wavelength (nm)')
    ylabel('Coefficient')
    legend({'Daylight PC1','Daylight PC2','Peak sensitivities'},'Location','best')
end

end