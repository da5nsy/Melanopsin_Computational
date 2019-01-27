clear, clc, close all;

try load('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\melcomp_1_results.mat')
catch
    
    %figure, hold on %to get figure to plot, also turn on 'plot_it' in melcomp
    
    range= [-230:10:-60,-59:135,140:10:390]; %94 seconds
    
    tic
    for i= 1:length(range)
        [MB1_minSD(i),MB2_minSD(i),melpeak(i),MB1_zeroSD(i),MB2_zeroSD(i),spread(:,i),MBx_m(:,:,i)]=melcomp_1(range(i));
        disp(melpeak(i))
        drawnow
    end
    toc
    
    save('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\results.mat')
end

%%

norm = 1;

figure, hold on


plot([melpeak(range==0),melpeak(range==0)],[0,max([MB1_zeroSD,MB2_zeroSD])],'k','LineWidth',4,'DisplayName','Mel Peak')
if norm
    plot(melpeak,ones(length(melpeak),1),'Color',[0.5,0.5,0.5],'LineWidth',4,'DisplayName','Baseline SD')
    scatter(melpeak,MB1_minSD/max(MB1_minSD),'r','filled','DisplayName','MB1: Min SD')
    scatter(melpeak,MB2_minSD/max(MB2_minSD),'b','filled','DisplayName','MB2: Min SD')
    
    
else
    plot(melpeak,MB1_zeroSD,'r','LineWidth',4,'DisplayName','MB1: Baseline SD')
    plot(melpeak,MB2_zeroSD,'b','LineWidth',4,'DisplayName','MB2: Baseline SD')
    scatter(melpeak,MB1_minSD,'r','filled','DisplayName','MB1: Min SD')
    scatter(melpeak,MB2_minSD,'b','filled','DisplayName','MB2: Min SD')
    %scatter(melpeak,spread(1,:),'r')
    %scatter(melpeak,spread(2,:),'b')
end

legend('Location','best')

xlabel('Hypothetical Melanopic Peak (nm)')
ylabel('Minimum standard deviation')

%% Plot one over the other

figure, hold on
scatter(melpeak,MB1_minSD/max(MB1_minSD),'r','filled')
scatter(melpeak,MB2_minSD/max(MB2_minSD),'b','filled')
title('Min SD for each object')

figure, hold on
scatter(melpeak,spread(1,:)/max(spread(1,:)),'r','filled')
scatter(melpeak,spread(2,:)/max(spread(2,:)),'b','filled')
title('SD of the set of means of objects')

figure, hold on
scatter(melpeak,spread(1,:)./MB1_minSD,'r','filled')
scatter(melpeak,spread(2,:)./MB2_minSD,'b','filled')
title('SD of set / Min SD for each object')

%% Plot the averages over wavelength
figure, hold on
xlim([0 2])
ylim([-4 4])
for i=1:size(MBx_m,3)-40
    scatter(MBx_m(1,:,i),MBx_m(2,:,i),'filled')
    %pause(0.1)
    drawnow
    title(melpeak(i))
end
