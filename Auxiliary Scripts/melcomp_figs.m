function melcomp_figs(plot_where,plot_size,base,ff,p)

% figures for paper

try %if within a function, do nothing
    nargin;
    if nargin ~= 3
        error
    end
catch %else, use default values
    clear, clc, close all
    
    plot_where = [20,60];
    plot_size  = [900,400];
    
    base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\figs';
    %base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\VNS18 Poster';
    
    ff = '-depsc'; %file format
    %ff = '-dpdf';
    
    p = 0; % 0 = figures display but don't save. 1 = figures display and save.
    
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

%% Fig: Mel_Phot

Mel_Phot_Corr

set(gcf,'Position',[plot_where plot_size],'color','w');

if p
    print([base,'\','Mel_Phot'],ff)
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
xticks(0:0.2:1); yticks(0:0.2:1);


subplot(1,2,2), axis square, hold on

plot([0.6,0.6],[0,0.93],'r')
plot([0,0.6],[0.033,0.033],'r:')
plot([0,0.6],[0.42,0.42],'r:')
plot([0,0.6],[0.93,0.93],'r:')

plot(y2,x,'k')
xlabel('Input'),ylabel('Output')
xticks(0:0.2:1); yticks(0:0.2:1);

if p
    print([base,'\','Monotonicity_concept'],ff)
end

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

if p
    print([base,'\','True3D'],ff)
end

%% Fig: chromaticities

figure('Position',[plot_where plot_size])
melcomp_2(1,1,1,1,'3D'); %Z-ax is arbitrary here
view(0,90)
xlim([0.5 1])

if p
    print([base,'\','chromaticities'],ff)
end

%% Fig: ZL

% Not currently used

figure('Position',[plot_where plot_size])
melcomp_2(1,1,1,9,'3D'); %final number sets Z-axis selection, 1 = 'L'
xlim([0.5 1])
view(340,40)
set(gcf,'color','w');

if p
    print([base,'\','ZL'],ff)
end

%% Fig: res_LMSRI

figure('Position',[plot_where plot_size.*[1,2.5]]) %bigger plot than standard
for i=1:5
    subplot(5,2,i*2-1)
    melcomp_2(1,1,1,i,'3D'); %final number sets Z-axis selection, 1 = 'L'
    view(0,0)
    xlim([0.6 0.8])
    if i ~= 5
        xlabel([])
        xticklabels([])
    end
end
for i=1:5
    subplot(5,2,i*2)
    melcomp_2(1,1,1,i,'3D'); %final number sets Z-axis selection, 1 = 'L'
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

if p
    print([base,'\','res_LMSRI'],ff)
end

%% Fig: res_lsri

figure('Position',[plot_where plot_size.*[1,2]]) %bigger plot than standard

for i=6:9
    subplot(4,2,(i-5)*2-1)
    melcomp_2(1,1,1,i,'3D'); %final number sets Z-axis selection, 1 = 'L'
    view(0,0)
    xlim([0.6 0.8])
    if i ~= 9
        xlabel([])
        xticklabels([])
    end
end
for i=6:9
    subplot(4,2,(i-5)*2)
    melcomp_2(1,1,1,i,'3D'); %final number sets Z-axis selection, 1 = 'L'
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

if p
    print([base,'\','res_lsri'],ff)
end

%% Fig: CTR

figure('Position',[plot_where plot_size.*[1,3]])
subplot(3,2,[1,4])
melcomp_2(1,1,1,9,'CTR');
text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,2,5)
melcomp_2(1,1,1,9,'CTR');
view(0,0)
legend('off')
text(0.02,0.98,'B','Units', 'Normalized', 'VerticalAlignment', 'Top')

subplot(3,2,6)
melcomp_2(1,1,1,9,'CTR');
view(90,0)
legend('off')
text(0.02,0.98,'C','Units', 'Normalized', 'VerticalAlignment', 'Top')

if p
    print([base,'\','CTR'],ff)
end

%% Fig: why_correl

melcomp_why('correl')

set(gcf,'Position',[plot_where plot_size]);

if p
    print([base,'\','why_correl'],ff)
end

%% Fig: why_PCA

melcomp_why('PCA')

set(gcf,'Position',[plot_where plot_size]);

if p
    print([base,'\','why_PCA'],ff)
end

%% Fig: inputs

%copied out of melcomp_2.m

figure('Position',[plot_where plot_size.*[1,2]])

subplot(3,1,1), hold on
D_CCT=1./linspace(1/3600,1/25000,20); %non-linear range, aiming to better reproduce observed variation
load B_cieday
T_SPD = GenerateCIEDay(D_CCT,[B_cieday]); %these appear to be linearly upsampled from 10nm intervals (see 'cieday investigation.m' https://github.com/da5nsy/General-Purpose-Functions/blob/3ee587429e9c4f3dd52d64acd95acf82d7e05f47/cieday%20investigation.m)
T_SPD = (T_SPD./repmat(max(T_SPD),81,1)); %normalise
S_SPD = S_cieday;
plot(SToWls(S_SPD),T_SPD)
xlim([350 800]); ylim([0 1]);
ylabel({'20 D-series illuminants'; 'Normalised SPD'});

subplot(3,1,2), hold on
load T_cones_sp
load T_rods
load T_melanopsin
plot(SToWls(S_cones_sp),T_cones_sp(3,:))
plot(SToWls(S_melanopsin),T_melanopsin)
plot(SToWls(S_rods),T_rods)
plot(SToWls(S_cones_sp),T_cones_sp(2,:))
plot(SToWls(S_cones_sp),T_cones_sp(1,:))

xlim([350 800]); ylim([0 1]);
ylabel({'5 opsins'; 'Normalised Spectral Sensitivity'});
legend('s-cone','melanopsin','rhodopsin','m-cone','l-cone');

subplot(3,1,3), hold on
load sur_vrhel
refs=[87, 93, 94, 134, 137, 138, 65, 19, 24, 140, 141];
T_refs = sur_vrhel';
T_refs = sur_vrhel(:,refs)';
S_refs = S_vrhel;
clear sur_vrhel refs S_vrhel
plot(SToWls(S_refs),T_refs)
xlim([350 800]); ylim([0 1]);
ylabel({'11 natural surfaces'; 'Normalised Spectral Reflectance'});

xlabel('Wavelength (nm)')

if p
    print([base,'\','inputs'],ff)
end

%% Fig: optimality
ss      = 1; %spectral sensitvities
prog    = 1; %progess read-out
t_range = [-200,2,400]; % target range, [start,interval,end]

figure('units','normalized','outerposition',[0 0 1 1]), hold on
%figure('Position',[plot_where plot_size]), hold on
melcomp_optimality(ss,prog,t_range)

%% Tight subplot demo

% clear, clc
%
% [ha, pos] = tight_subplot(3,2);
% for ii = 1:6;
%     axes(ha(ii));
%     plot(randn(10,ii));
% end
% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')