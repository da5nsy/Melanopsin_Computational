clear, clc, close all

load melcomp_3_fullWorkspace.mat
load melcomp_3_correlation_results.mat

cs_b = cs(:,:,squeeze(mean(cs(1,:,:),2))>0.5); %cs_bright. Same criterea as curve fitting. Arbitrary threshold. Here only applies to L, not ideal.


%% Plot visualisation of how close estimates of PC2 are from different signals

plt_pc2estimates = 0;

if plt_pc2estimates
    reflectances = [1:5:size(cs_b,2)];
    illuminants = [1:100:size(cs_b,3)];
    
    figure, hold on
    
    for illum = 150%illuminants
        
        plot3([0,max(reflectances)],[pc_p.score(illum,2),pc_p.score(illum,2)],[illum illum],'k:')
        
        for ref = 1:size(cs_b,2)%reflectances
            
            fv = 14;
            sv = 3;
            
            y = (p_1(fv,sv,1) * cs_b(sv,ref,illum) + p_1(fv,sv,2))...
                * cs_b(fv,ref,illum) + ...
                (p_2(fv,sv,1) * cs_b(sv,ref,illum) + p_2(fv,sv,2));
            scatter3(ref,y,illum,'g')
            
            fv = 18;
            sv = 1;
            
            y = (p_1(fv,sv,1) * cs_b(sv,ref,illum) + p_1(fv,sv,2))...
                * cs_b(fv,ref,illum) + ...
                (p_2(fv,sv,1) * cs_b(sv,ref,illum) + p_2(fv,sv,2));
            scatter3(ref,y,illum,'r')
        end
    end
    
    xlabel('Reflectance Sample')
    ylabel('Predicted PC2')
    zlabel('Illuminant')
    legend({'line','14,3','18,1'})
    
    xticks(1:120)
    
    load sur_vrhel_withLabels.mat
    %labels_vrhel_nat = [labels_vrhel([1:44,65,69,81:154]).label];
    xticklabels([labels_vrhel([1:44,65,69,81:154]).label])
    set(gca,'XTickLabelRotation',90)
end



%%

for i=1:21
    for j=1:length(cs_b)
        a(i,j) = std(cs_b(i,:,j));
    end
end

%%
figure, hold on
for i=1:21
    plot3(1:100,a(i,1:100),ones(100)*i)
end

xlabel('Illumination #')
ylabel('standard deviation of signals')
zlabel('signal')

ylim([0 0.2])