clear, clc, close all

plot_where = [500,200];
plot_size  = [800,400];

set(groot,'defaultfigureposition',[100 100 500 400]); 
set(groot,'defaultLineLineWidth',2);
set(groot,'defaultAxesFontName', 'Courier');
set(groot,'defaultAxesFontSize',12);
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vectors
set(groot,'defaultfigurecolor','white')

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_1_caller';

print_figures = 1;

%%

try load('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\melcomp_1_results.mat')
catch
    
    figure, hold on
    
    range= [-150:10:-60,-58:2:138,140:10:200];
    
    MB1_minSD = zeros(1,length(range));
    MB2_minSD = zeros(1,length(range));
    spread =    zeros(2,length(range));
    MBx_m  =    zeros(2,11,length(range));
    
    tic
    for i= 1:length(range)
        [MB1_minSD(i),MB2_minSD(i),MB1_zeroSD,MB2_zeroSD,spread(:,i),MBx_m(:,:,i)]=melcomp_1(range(i),'plt_it_overide',1);
        disp(range(i))
        drawnow
    end
    toc
    
    save('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\melcomp_1_results.mat')
    
    title([]);
    ylim([0 max([MB1_zeroSD,MB2_zeroSD])])
    xticks([min(xlim),0,max(xlim)])
    yticks([0, max(ylim)])
    
    if print_figures
        save2pdf([base,'\it.pdf'])
        savefig([base,'\it.fig'])
    end
end

melpeak = 488+range;

%% Minimals SD

norm = 1;

figure, hold on

if norm
    plot([melpeak(range==0),melpeak(range==0)],[0,1],...
        'k','LineWidth',4,'DisplayName','Mel Peak')
    plot(melpeak,ones(length(melpeak),1),'Color',[0.5,0.5,0.5],'LineWidth',4,'DisplayName','Baseline SD')
    scatter(melpeak,MB1_minSD/max(MB1_minSD),'r','filled','DisplayName','MB1')
    scatter(melpeak,MB2_minSD/max(MB2_minSD),'b','filled','DisplayName','MB2')
    
    
else
    plot([melpeak(range==0),melpeak(range==0)],[0,max([MB1_zeroSD,MB2_zeroSD])],...
        'k','LineWidth',4,'DisplayName','Mel Peak')
    plot([melpeak(1),melpeak(end)],[MB1_zeroSD,MB1_zeroSD],'r','LineWidth',4,'DisplayName','MB1: Baseline SD')
    plot([melpeak(1),melpeak(end)],[MB2_zeroSD,MB2_zeroSD],'b','LineWidth',4,'DisplayName','MB2: Baseline SD')
    scatter(melpeak,MB1_minSD,'r','filled','DisplayName','MB1: Min SD')
    scatter(melpeak,MB2_minSD,'b','filled','DisplayName','MB2: Min SD')
    %scatter(melpeak,spread(1,:),'r')
    %scatter(melpeak,spread(2,:),'b')
end

legend('Location','best')
xlabel('Nominal Melanopic Peak (nm)')
ylabel({'Minimum standard deviation',' of whole set (normalised)'})
axis tight
yticks([0 1])

if print_figures
    save2pdf([base,'\opt.pdf'])
end

%% Spread

figure, hold on
scatter(melpeak,spread(1,:)/max(spread(1,:)),'r','filled')
scatter(melpeak,spread(2,:)/max(spread(2,:)),'b','filled')
xlabel('Nominal Melanopic Peak (nm)')
ylabel('Mean of inter-object distances')
axis tight
ylim([0 1])
yticks([0 1])
if print_figures
    save2pdf([base,'\sdmeans.pdf'])
end

%% Plot the averages over wavelength

figure('Position',[plot_where plot_size]), hold on
xlim([0.3 1.1])
ylim([-0.05 0.15])
col2 = parula(size(MBx_m,3));

for i = 1:size(MBx_m,3)
   chull = convhull(MBx_m(1,:,i),MBx_m(2,:,i));
   fill(MBx_m(1,chull,i),MBx_m(2,chull,i),col2(i,:),'FaceAlpha',0.6,'LineStyle','none')
end

colormap('parula')
colorbar
caxis([min(melpeak) max(melpeak)])

if print_figures
    save2pdf([base,'\space.pdf'])
end

%% - %% Visualiser

clear, clc, close all;

for i=-100:10:100
    melcomp_1(i,'plt_appf_overide',1);
    axis auto
    axis equal
    title(i)
    drawnow
    %pause(0.5)
end





















