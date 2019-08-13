% tests whether the spectral sensitivity of melanopsin is particularly good
% for this or not compared to other SSFs

range = -100:1:200; %with 488nm as the origin

for rep = 1:10
    for rI = 1:length(range) %range index
        KMM(rI,:,rep) = optimality(range(rI));
        disp(range(rI))
    end
end

%% Plot

figure, 
hold on

ylim([0 1])
%plot(488+range,squeeze(KMM(:,4,:)))
plot([488,488],ylim,'k','DisplayName','melanopsin')
plot(488+range,min(squeeze(KMM(:,4,:)),[],2),'r','DisplayName','min')
plot(488+range,mean(squeeze(KMM(:,4,:)),2),'g','DisplayName','mean')
plot(488+range,max(squeeze(KMM(:,4,:)),[],2),'b','DisplayName','max')
legend('Location','best')

xlabel('Nominal mel peak SSF (nm)')
ylabel('k-means-mark')

%% Save fig

save2pdf('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\optimality_caller\optimality.pdf')
