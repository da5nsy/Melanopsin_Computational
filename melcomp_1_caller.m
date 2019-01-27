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
    scatter(melpeak,MB1_minSD/max(MB1_minSD),'r','filled','DisplayName','MB1')
    scatter(melpeak,MB2_minSD/max(MB2_minSD),'b','filled','DisplayName','MB2')
    
    
else
    plot(melpeak,MB1_zeroSD,'r','LineWidth',4,'DisplayName','MB1: Baseline SD')
    plot(melpeak,MB2_zeroSD,'b','LineWidth',4,'DisplayName','MB2: Baseline SD')
    scatter(melpeak,MB1_minSD,'r','filled','DisplayName','MB1: Min SD')
    scatter(melpeak,MB2_minSD,'b','filled','DisplayName','MB2: Min SD')
    %scatter(melpeak,spread(1,:),'r')
    %scatter(melpeak,spread(2,:),'b')
end

legend('Location','best')

xlabel('Nominal Melanopic Peak (nm)')
ylabel('Minimum standard deviation of whole set (normalised)')

%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\Melanopsin Computational\figs\opt.pdf")

%% Plot one over the other

figure, hold on
scatter(melpeak,spread(1,:)/max(spread(1,:)),'r','filled','DisplayName','MB1')
scatter(melpeak,spread(2,:)/max(spread(2,:)),'b','filled','DisplayName','MB2')
legend('Location','best')
xlabel('Nominal Melanopic Peak (nm)')
ylabel('SD of the set of means of objects (normalised)')
%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\Melanopsin Computational\figs\sdmeans.pdf")

figure, hold on
scatter(melpeak,(spread(1,:)/max(spread(1,:)))./(MB1_minSD/max(MB1_minSD)),'r','filled','DisplayName','MB1')
scatter(melpeak,(spread(2,:)/max(spread(2,:)))./(MB2_minSD/max(MB2_minSD)),'b','filled','DisplayName','MB2')
legend('Location','best')
xlabel('Nominal Melanopic Peak (nm)')
ylabel('SD of set / Min SD for each object')
%save2pdf("C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Thesis\Melanopsin Computational\figs\setoverobj.pdf")

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

%% - %% Visualiser

clear, clc, %close all;

for i=-100:10:100
    melcomp_1(i)
    axis auto
    axis equal
    title(i)
    drawnow
    %pause(0.5)
end



