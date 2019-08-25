clear, clc, close all

%Pre-flight
d.s=25;               % display, size
d.MFA = 0.2;          % Marker Face Alpha
d.mktrns = 0.3;       % Marker transparency
set(groot,'defaultfigureposition',[100 100 500 400])
set(groot,'defaultLineLineWidth',2)
set(groot,'defaultAxesFontName', 'Courier')
set(groot,'defaultAxesFontSize',12)
set(groot,'defaultFigureRenderer', 'painters') %renders pdfs as vector graphics
set(groot,'defaultfigurecolor','white')
cols = hsv(10); rng(2);
set(groot,'defaultAxesColorOrder',cols(randperm(size(cols,1)),:))

base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\KMeansMarkDemo';

plt.print = 0;
if plt.print
    warning('plt.print is enabled - you sure? This will overwrite existing figures.')
end


%%
wholeset = 1;
[~,~,~,~,MB_star1] = transformToIllIndSpace(0,wholeset);
[~,~,~,~,MB_star2] = transformToIllIndSpace(590-488,wholeset);

close all

%%

[KMM1,km_idx1] = KMeansMark(MB_star1);
disp(KMM1)
figure,
scatter(MB_star1(1,:),MB_star1(2,:),d.s,cols(km_idx1,:),'filled','MarkerFaceAlpha',d.MFA)
legend off
cleanTicks
xlabel('{\itl}_{MB} + {\itk_1i}_{MB}');
ylabel('{\its}_{MB} + {\itk_2i}_{MB}');

if plt.print
    save2pdf([base,'\1.pdf'])
end

[KMM2,km_idx2] = KMeansMark(MB_star2);
disp(KMM2)
figure,
scatter(MB_star2(1,:),MB_star2(2,:),d.s,cols(km_idx2,:),'filled','MarkerFaceAlpha',d.MFA)
legend off
cleanTicks
xlabel('{\itl}_{MB} + {\itk_1i}_{MB}');
ylabel('{\its}_{MB} + {\itk_2i}_{MB}');

if plt.print
    save2pdf([base,'\2.pdf'])
end

%%
clc % lots of warning messages from kmeans
figure, hold on
for n = 1:15
    rng(4)
    out(:,n) = kmeans(MB_star1(:,:)',10,'MaxIter',n);
    cla
    %gscatter(MB_star1(1,:),MB_star1(2,:),out(:,n))    
    scatter(MB_star1(1,:),MB_star1(2,:),d.s,cols(out(:,n),:),'filled','MarkerFaceAlpha',d.MFA)
    for i = 1:10
        scatter(mean(MB_star1(1,out(:,n) == i)),mean(MB_star1(2,out(:,n) == i)),'k*')
    end
    xlabel('{\itl}_{MB} + {\itk_1i}_{MB}');
    ylabel('{\its}_{MB} + {\itk_2i}_{MB}');
    title(n)
    drawnow      
    pause(0.1)
end


