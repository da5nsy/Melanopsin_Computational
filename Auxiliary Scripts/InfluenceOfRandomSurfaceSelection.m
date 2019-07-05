% Running multiple versions of melcomp_9 to understand the variability

clear, clc, close all

for i=1:10
    mark(:,:,i) = melcomp_9;
    close all
    disp(['RUN ',num2str(i)])
end

%%

titles = {'Do nothing','Grey World','Bright-is-White','Melanopsin'};
markers = {'s:','d:','^:','o:'};
mfc = hsv(4);
lw = 3;

figure, hold on
for j=1:10
for i=1:4
    plot([12:7:100,100],mark(i,:,j),markers{i},'Color',mfc(i,:),'MarkerFaceColor',mfc(i,:),'linewidth',lw)
end
end
legend(titles,'Location','Northwest')
%plot(p.Results.pcSurfRange,1./round((p.Results.pcSurfRange/100) * size(T_SRF_reduced,2)),'k','DisplayName','Chance','linewidth',lw-1)
xlim([0 100])
ylim([0 1])

xlabel(sprintf('Percentage of surfaces used in each run (/%d)',size(T_SRF_reduced,2)))
ylabel('K-means-mark (/1)')
grid on