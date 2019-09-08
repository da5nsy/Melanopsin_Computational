clear, clc, close all
d = DGdisplaydefaults;

[T_SPD, ~, T_SSF, ~, S_sh] = melcomp_loader('SRF','Vrhel_nat_extended','SSF','SP','Lum','SP');

%% - %%

% Copied from melomp_3 :

%% PCA of daylight

% PCA of spectral Power distributions
% compute pca variable weight
vw = 2;
if vw == 0 %no variable weights
    pc_p.variableweights = ones(81,1);
elseif vw == 1 %S + L and max in between
    pc_p.variableweights = T_SSF(:,1)+T_SSF(:,3);
    [~,t_p_locs] = findpeaks(pc_p.variableweights);
    pc_p.variableweights(t_p_locs(1):t_p_locs(2)) = max(pc_p.variableweights);
    pc_p.variableweights(pc_p.variableweights>1) = 1;
elseif vw == 2 %S+M+L+I
    pc_p.variableweights = sum(T_SSF');
end

plt_vw = 1;
if plt_vw
    figure,
    axis tight
    plot(SToWls(S_sh),pc_p.variableweights)
    xlim([390 730]);
    xticks(400:100:700);
    %ylim([0 1]);
    yticks(ylim);
    ylabel('Weight')
end

[pc_p.coeff, pc_p.score, pc_p.latent, pc_p.tsquared, pc_p.explained, pc_p.mu] = pca(T_SPD','VariableWeights',pc_p.variableweights);

plt_pc_p = 1;
if plt_pc_p
    figure, hold on
    axis tight
    xlim([390 730]); xticks(400:100:700);
    ylim([-0.55 0.55]); yticks([-0.5,0,0.5]);
    %legend({'PC1','PC2','PC3'},'Location','Southwest')
    
    % highlights negative space
    fill([min(xlim),max(xlim),max(xlim),min(xlim),min(xlim)],[-1,-1,0,0,-1],...
        'k','LineStyle','none','FaceAlpha','0.03');
    
    %normalised PCs:
    %plot(SToWls(S_sh),pc.coeff(:,1:3)./max(pc.coeff(:,1:3)))
    
    plot(SToWls(S_sh),pc_p.coeff(:,1:3))    
   
end



