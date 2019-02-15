function melcomp_2_figs(plot_where,plot_size,base,p)

% figures for paper

try %if within a function, do nothing
    nargin;
    if nargin ~= 3
        error
    end
catch %else, use default values
    clear, clc, close all
    
    plot_where = [500,200];
    plot_size  = [800,400];
    
    %base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\figs';
    base = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\melcomp_2_figs';

    p = 1; % 0 = figures display but don't save. 1 = figures display and save.
    set(0,'defaultAxesFontName', 'Courier')
    
    disp('using default values') %!!!!!!!!!!!!!! this is not displaying currently !!!!!!!!!!!!!!!!
end

% melcomp_2.m reference notes:
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

% expectedSPD = {'Granada_sub','Granada','D-series'};
% expectedSRF = {'Vrhel_nat_1','Vrhel_nat_2','Vrhel_full','Foster'};
% expectedSSF = {'SS10','SP'};
% expectedlum = {'CIE_10','SP'};
% default_SPD = expectedSPD{3};
% default_SRF = expectedSRF{1};
% default_SSF = expectedSSF{2};
% default_lum = expectedlum{2};

%% Fig: Mel_Phot

Mel_Phot_Corr

set(gcf,'Position',[plot_where plot_size],'color','w');
cleanTicks

if p
    save2pdf([base,'\Mel_Phot.pdf'])
end

%% Fig: Monotonicity_concept

x  = 0:0.01:1;
y1 = (x.^2)/max(x.^2);
y2 = (sin(x*7)+1)/2;

figure('Position',[plot_where plot_size])

subplot(1,2,1), axis square, hold on

plot([0.6,0.6],[0,y1(61)],'r')
plot([0,0.6],[y1(61),y1(61)],'r')

plot(x,y1,'k')

xlabel('Input'),ylabel('Output')
cleanTicks


subplot(1,2,2), axis square, hold on

plot([0.6,0.6],[0,0.93],'r')
plot([0,0.6],[0.033,0.033],'r:')
plot([0,0.6],[0.42,0.42],'r:')
plot([0,0.6],[0.93,0.93],'r:')

plot(y2,x,'k')
xlabel('Input'),ylabel('Output')
cleanTicks

if p
    save2pdf([base,'\Monotonicity_concept.pdf'])
end

%% Fig: True3D

rng(4) %Generate same 'random' numbers each time
x=rand(10,1);
y=rand(10,1);

figure('Position',[plot_where plot_size])

subplot(1,2,1), axis square
scatter(x,y,'k','filled');
xlim([0 1]); ylim([0 1]);
axis square
cleanTicks
xlabel('x'); ylabel('y');

subplot(1,2,2), axis square
scatter(x,y./y,'k','filled');
xlim([0 1]); ylim([0 1]);
axis square
cleanTicks
xlabel('x'); ylabel('y/y');

if p
    save2pdf([base,'\True3D.pdf'])
end

%% Fig: chromaticities

figure('Position',[plot_where plot_size])
melcomp_2('plt','3D');
view(0,90)
xlim([0.5 1])
cleanTicks

if p
    save2pdf([base,'\chromaticities.pdf'])
end

%% Fig: ZL

% Not currently used

figure('Position',[plot_where plot_size])
melcomp_2('plt','3D'); %final number sets Z-axis selection, 1 = 'L'
xlim([0.5 1])
view(340,40)
set(gcf,'color','w');
cleanTicks

if p
    save2pdf([base,'\ZL.pdf'])
end

%% Fig: res_LMSRI

figure('Position',[plot_where plot_size.*[1,2.5]]) %bigger plot than standard
for i=1:5
    subplot(5,2,i*2-1)
    melcomp_2('Z_ax',i,'plt','3D'); 
    view(0,0)
    xlim([0.6 0.8])
    cleanTicks
    if i ~= 5
        xlabel([])
        xticklabels([])
    end
end
for i=1:5
    subplot(5,2,i*2)
    melcomp_2('Z_ax',i,'plt','3D');
    view(90,0)
    ylim([0 0.025])
    cleanTicks
    if i ~= 5
        ylabel([])
        yticklabels([])
    end
    zticklabels([])
    zlabel([])
end

set(gcf,'color','w');

if p    
    save2pdf([base,'\res_LMSRI.pdf'])
end

%% Fig: res_lsri

figure('Position',[plot_where plot_size.*[1,2]]) %bigger plot than standard

for i=6:9
    subplot(4,2,(i-5)*2-1)
    melcomp_2('Z_ax',i,'plt','3D'); %final number sets Z-axis selection, 1 = 'L'
    view(0,0)
    xlim([0.6 0.8])
    cleanTicks
    if i ~= 9
        xlabel([])
        xticklabels([])
    end
end
for i=6:9
    subplot(4,2,(i-5)*2)
    melcomp_2('Z_ax',i,'plt','3D'); %final number sets Z-axis selection, 1 = 'L'
    view(90,0)
    ylim([0 0.025])
    cleanTicks
    if i ~= 9
        ylabel([])
        yticklabels([])
    end
    zticklabels([])
    zlabel([])
end

set(gcf,'color','w');

if p        
    save2pdf([base,'\res_lsri.pdf'])
end

%% Fig: CTR

figure('Position',[plot_where plot_size.*[1,3]])
subplot(3,2,[1,4])
melcomp_2('plt','CTR');
text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top')
cleanTicks

subplot(3,2,5)
melcomp_2('plt','CTR');
view(0,0)
legend('off')
text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top')
cleanTicks

subplot(3,2,6)
melcomp_2('plt','CTR');
view(90,0)
legend('off')
text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top')
cleanTicks

if p        
    save2pdf([base,'\CTR.pdf'])
end

%% Fig: why_correl

melcomp_why('correl')

set(gcf,'Position',[plot_where plot_size]);

if p
    save2pdf([base,'\why_correl.pdf'])
end

%% Fig: why_PCA

melcomp_why('PCA')

set(gcf,'Position',[plot_where plot_size]);

if p    
    save2pdf([base,'\why_PCA.pdf'])
end

%% Fig: inputs

%copied out of melcomp_2.m

figure('Position',[plot_where plot_size.*[1,2]])

[T_SPD, T_SRF, T_SSF, ~, S_sh] = melcomp_loader('SPD','D-series','SRF','Vrhel_nat_1','SSF','SP','lum','SP');


subplot(3,1,1), hold on
plot(SToWls(S_sh),T_SPD)
xlim([S_sh(1) S_sh(1)+S_sh(2)*(S_sh(3)-1)]); ylim([0 1]);
ylabel({'20 D-series illuminants'; 'Normalised SPD'});
cleanTicks

subplot(3,1,2), hold on
plot(SToWls(S_sh),T_SSF(:,3))
plot(SToWls(S_sh),T_SSF(:,5))
plot(SToWls(S_sh),T_SSF(:,4))
plot(SToWls(S_sh),T_SSF(:,2))
plot(SToWls(S_sh),T_SSF(:,1))
xlim([S_sh(1) S_sh(1)+S_sh(2)*(S_sh(3)-1)]); ylim([0 1]);
ylabel({'5 opsins'; 'Normalised Spectral Sensitivity'});
legend('s-cone','melanopsin','rhodopsin','m-cone','l-cone');
cleanTicks

subplot(3,1,3), hold on
plot(SToWls(S_sh),T_SRF)
xlim([S_sh(1) S_sh(1)+S_sh(2)*(S_sh(3)-1)]); ylim([0 1]);
ylabel({'11 natural surfaces'; 'Normalised Spectral Reflectance'});
cleanTicks

xlabel('Wavelength (nm)')

if p    
    save2pdf([base,'\inputs.pdf'])
end

%% Fig: optimality 
% not currently used
% 
% ss      = 1; %spectral sensitvities
% prog    = 1; %progess read-out
% t_range = [-180,2,400]; % target range, [start,interval,end]
% 
% figure('units','normalized','outerposition',[0 0 1 1]), hold on
% %figure('Position',[plot_where plot_size]), hold on
% melcomp_optimality(ss,prog,t_range)
