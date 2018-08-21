% figures for paper

clc, clear, close all

plot_where = [20,60];
plot_size  = [900,400];

% melcomp.m reference notes:
% PF_SPD = 1;
% % 1 = CIE D series
% % 2 = Hernández-Andrés+
% 
% PF_refs = 1;
% % 1 = Vhrel+ (natural only)
% % 2 = Vhrel+ (all)
% % 3 = Foster+
% 
% PF_obs = 1;
% % 1 = PTB Smith-Pokorny

%

% plt_lbls{1}  = 'L'; 
% plt_lbls{2}  = 'M';
% plt_lbls{3}  = 'S';
% plt_lbls{4}  = 'R';
% plt_lbls{5}  = 'I';
% plt_lbls{6}  = 'l';
% plt_lbls{7}  = 's';
% plt_lbls{8}  = 'r';
% plt_lbls{9}  = 'i';
% plt_lbls{10} = 'L+M';
% plt_lbls{11} = '(0.6373*L)+(0.3924*M)';
% plt_lbls{12} = 'r + i';

%% Fig: Mel_Phot

Mel_Phot_Corr

set(gcf,'Position',[plot_where plot_size],'color','w');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\Mel_Phot','-depsc')

%% Fig: Monotonicity_concept

x  = 0:0.01:1;
y1 = (x.^2)/max(x.^2);
y2 = (sin(x*7)+1)/2;

figure('Position',[plot_where plot_size])

subplot(1,2,1), axis square, hold on

plot([0.6,0.6],[0,y1(61)],'r')
plot([0,0.6],[y1(61),y1(61)],'r')

plot(x,y1,'k')

xlabel('Input'),ylabel('Determined Output')
xticks(0:0.2:1); yticks(0:0.2:1);


subplot(1,2,2), axis square, hold on

plot([0.6,0.6],[0,0.93],'r')
plot([0,0.6],[0.033,0.033],'r:')
plot([0,0.6],[0.42,0.42],'r:')
plot([0,0.6],[0.93,0.93],'r:')

plot(y2,x,'k')
xlabel('Input'),ylabel('Undetermined Output')
xticks(0:0.2:1); yticks(0:0.2:1);

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\Monotonicity_concept','-depsc')

%% Fig: True3D

rng(4) %Generate same 'random' numbers each time
x=rand(10,1);
y=rand(10,1);

figure('Position',[plot_where plot_size])

subplot(1,2,1), axis square
scatter(x,y,'k','filled');
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('x'); ylabel('y');

subplot(1,2,2), axis square
scatter(x,y./y,'k','filled');
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('x'); ylabel('y/y');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\True3D','-depsc')

%% Fig: chromaticities

figure('Position',[plot_where plot_size])
melcomp(1,1,1,1,'3D') %Z-ax is arbitrary here
view(0,90)
xlim([0.5 1])

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\chromaticities','-depsc')


%% Fig: ZL

% Not currently used

figure('Position',[plot_where plot_size])
melcomp(1,1,1,1,'3D') %final number sets Z-axis selection, 1 = 'L'
xlim([0.5 1])
view(340,40)
set(gcf,'color','w');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\ZL','-depsc')

%% Fig: res_LMSRI

figure('Position',[plot_where plot_size.*[1,2.5]]) %bigger plot than standard
for i=1:5
    subplot(5,2,i*2-1)
    melcomp(1,1,1,i,'3D') %final number sets Z-axis selection, 1 = 'L'
    view(0,0)    
    xlim([0.6 0.8])
    if i ~= 5
        xlabel([])
        xticklabels([])
    end
end
for i=1:5
    subplot(5,2,i*2)
    melcomp(1,1,1,i,'3D') %final number sets Z-axis selection, 1 = 'L'
    view(90,0)
    ylim([0 0.8])
    if i ~= 5
        ylabel([])
        yticklabels([])
    end
    zticklabels([])
    zlabel([])
end

set(gcf,'color','w');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\res_LMSRI','-depsc')

%% Fig: res_lsri

figure('Position',[plot_where plot_size.*[1,2]]) %bigger plot than standard

for i=6:9
    subplot(4,2,(i-5)*2-1)
    melcomp(1,1,1,i,'3D') %final number sets Z-axis selection, 1 = 'L'
    view(0,0)    
    xlim([0.6 0.8])
    if i ~= 9
        xlabel([])
        xticklabels([])
    end
end
for i=6:9
    subplot(4,2,(i-5)*2)
    melcomp(1,1,1,i,'3D') %final number sets Z-axis selection, 1 = 'L'
    view(90,0)
    ylim([0 0.8])
    if i ~= 9
        ylabel([])
        yticklabels([])
    end
    zticklabels([])
    zlabel([])
end

set(gcf,'color','w');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\res_lsri','-depsc')

%% Fig: CTR

figure('Position',[plot_where plot_size.*[1,3]])
subplot(3,2,[1,4])
melcomp(1,1,1,9,'CTR')
text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,2,5)
melcomp(1,1,1,9,'CTR')
view(0,0)
legend('off')
text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,2,6)
melcomp(1,1,1,9,'CTR')
view(90,0)
legend('off')
text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top')

print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\CTR','-depsc')


%% Fig: why_correl

melcomp_why('correl')

set(gcf,'Position',[plot_where plot_size]);
%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\why_correl','-depsc')

%% Fig: why_PCA

melcomp_why('PCA')

set(gcf,'Position',[plot_where plot_size]);
%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\why_PCA','-depsc')

%% Tight subplot demo

% clear, clc
% 
% [ha, pos] = tight_subplot(3,2);
% for ii = 1:6;
%     axes(ha(ii));
%     plot(randn(10,ii));
% end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')