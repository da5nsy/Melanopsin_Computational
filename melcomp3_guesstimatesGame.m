%% Guesstimates game

% clearvars -except cs p_1 p_2 pc_p
% save('melcomp_3_gamedata.mat','cs','p_1','p_2','pc_p')

clc, clear, close all
load melcomp_3_gamedata.mat

%%

MorC = 1;

for k=1:size(cs,3)
    
    % OK, here's the game:
    % You get the cs for n surfaces under 1 illuminant.
    % Tell me what the PC2 score for that illuminant is.
    
    g_nSur=5; % number of surfaces
    
    %rng(4);
    g_ref_ind = randi(size(cs,2),g_nSur,1); %surfaces index
    g_ill_ind = k;% randi(size(cs,2),1,1)
    
    g_cs = cs(1:5,g_ref_ind,g_ill_ind);
    
    for i=1:5
        for j=1:5
            if i == j
                continue
            end
            g_fits(:,i,j) = orthogonalRegress(log10(g_cs(i,:)),log10(g_cs(j,:)));
            
            
            m = p_1(1,i,j) * (pc_p.score(g_ill_ind,1)-min(pc_p.score(:,1))) + p_1(2,i,j);
            c = p_2(1,i,j) * (pc_p.score(g_ill_ind,1)-min(pc_p.score(:,1))) + p_2(2,i,j);
            
            g_estimatedPC2(i,j) =  m .* g_fits(MorC,i,j) + c;
        end
    end
    
    correct_answer = pc_p.score(g_ill_ind,2);
    
    %
    i=3;
    j=5;
    
    % figure, hold on
    % scatter(log10(g_cs(i,:)),log10(g_cs(j,:)))
    %
    % xlin = linspace(min(log10(g_cs(i,:))),max(log10(g_cs(i,:))));
    % ylin = g_fits(1,i,j)*xlin+g_fits(2,i,j);
    % plot(xlin,ylin)
    
    %
    
    g_estimatedPC2_diff(k,:,:) = abs(g_estimatedPC2 - correct_answer);
    % figure,
    % imagesc(g_estimatedPC2_diff)
    % colormap('gray')
    
end

%
a = squeeze(mean(g_estimatedPC2_diff,1))

figure,
imagesc(a)
colormap('gray')

colorbar

axis image
% xticks(1:5); xticklabels(plt_lbls{1:5});
% yticks(1:5); yticklabels(plt_lbls{1:5});
% 
