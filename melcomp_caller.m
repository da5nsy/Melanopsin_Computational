clear, clc, close all;

%figure, hold on %to get figure to plot, also turn on 'plot_it' in melcomp

range= [-230:10:-60,-59:135,140:10:390];

for i= 1:length(range)
    [MB1_minSD(i),MB2_minSD(i),melpeak(i),MB1_zeroSD(i),MB2_zeroSD(i)]=melcomp(range(i));
    disp(range(i))
    drawnow
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
    
    legend({'488nm nominal peak of melanopsin',...
        'MB1: Baseline SD (no CA model)',...
        'MB2: Baseline SD (no CA model)',...
        'MB1: Minimum achievable SD with correction',...
        'MB2: Minimum achievable SD with correction'},...
        'Location','best')
end

xlabel('Hypothetical Melanopic Peak (nm)')
ylabel('Minimum standard deviation')





%%