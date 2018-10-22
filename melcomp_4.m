clear, clc, close all

load melcomp_3_fullWorkspace.mat
load melcomp_3_correlation_results.mat

%%
n_smp = 20; %number of samples
light = 500;

figure, hold on

for light = 500:520
    
    plot([0,n_smp],[pc_p.score(light,2),pc_p.score(light,2)],'k:')
    
    for smp=1:n_smp
        
        fv = 7;
        sv = 1;
        
        y = (p_1(fv,sv,1) * cs(sv,smp,light) + p_1(fv,sv,2))...
            * cs(fv,smp,light) + ...
            (p_2(fv,sv,1) * cs(sv,smp,light) + p_2(fv,sv,2));
        scatter(smp,y,'g')
        
        fv = 18;
        sv = 1;
        
        y = (p_1(fv,sv,1) * cs(sv,smp,light) + p_1(fv,sv,2))...
            * cs(fv,smp,light) + ...
            (p_2(fv,sv,1) * cs(sv,smp,light) + p_2(fv,sv,2));
        scatter(smp,y,'r')
    end
end














%%

% for i=1:21
%     for j=1:2600
%         a(i,j) = std(cs(i,:,j));
%     end
% end
%
% %%
% figure, hold on
% for i=1:21
%     plot3(1:100,a(i,1:100),ones(100)*i)
% end
%
% xlabel('Illumination #')
% ylabel('standard deviation of signals')
% zlabel('signal')
%
% ylim([0 0.2])