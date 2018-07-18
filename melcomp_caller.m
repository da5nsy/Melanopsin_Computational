clear, clc, close all;

try load('results.mat')
catch

    %figure, hold on %to get figure to plot, also turn on 'plot_it' in melcomp

    range= [-230:10:-60,-59:135,140:10:390]; %94 seconds

    tic
    for i= 1:length(range)
        [MB1_minSD(i),MB2_minSD(i),melpeak(i),MB1_zeroSD(i),MB2_zeroSD(i),spread(:,i),MBx_m(:,:,i)]=melcomp(range(i));
        disp(melpeak(i))
        drawnow
    end
    toc

    save('results.mat')
end

%%

norm = 0;

figure, hold on


plot([melpeak(range==0),melpeak(range==0)],[0,max([MB1_zeroSD,MB2_zeroSD])],'k','LineWidth',4)
if norm
    plot(melpeak,ones(length(melpeak),1),'k','LineWidth',4)
    scatter(melpeak,MB1_minSD/max(MB1_minSD),'r','filled')
    scatter(melpeak,MB2_minSD/max(MB2_minSD),'b','filled')
    
    
else
    plot(melpeak,MB1_zeroSD,'r','LineWidth',4)
    plot(melpeak,MB2_zeroSD,'b','LineWidth',4)
    scatter(melpeak,MB1_minSD,'r','filled')
    scatter(melpeak,MB2_minSD,'b','filled')
    scatter(melpeak,spread(1,:),'r')
    scatter(melpeak,spread(2,:),'b')
    
    legend({'488nm nominal peak of melanopsin',...
        'MB1: Baseline SD (no CA model)',...
        'MB2: Baseline SD (no CA model)',...
        'MB1: Minimum achievable SD with correction',...
        'MB2: Minimum achievable SD with correction'},...
        'Location','best')
end

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
    