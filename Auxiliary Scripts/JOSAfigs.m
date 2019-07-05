clear, clc, close all

plot_where = [500,300];
plot_size  = [1000,350];

min_l_scale = 0.6;
max_l_scale = 0.8;
max_s_scale = 0.04;
mkrsz = 15; %marker size
mktrns = 0.3; %marker transparency

set(0,'defaultAxesFontName', 'Courier')

base = 'C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\JOSA\figs';
print_figures = 0;

%%

% Load data
[T_SPD, T_SRF, T_SSF, T_lum, S_sh] = melcomp_loader(...
    'SPD','Granada_sub',...
    'SRF','Vrhel_nat_1',...
    'SSF','SS10',...
    'mel_offset',0,...
    'lum','CIE_10');

% Compute colorimetry
[LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh);

%lsri = log(lsri);

%%

figure('Position',[plot_where plot_size]),

sp = tight_subplot(1,3,[.01 .03],[.18 .05],[.07 .015]);

axes(sp(1))
scatter(lsri(1,:,1),lsri(2,:,1),mkrsz,'k','MarkerEdgeAlpha',mktrns)
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
cleanTicks
xlabel('{\itl}_{MB,10}');
ylabel('{\its}_{MB,10}');
text(0.05,0.95,'A','Units','normalized')

axes(sp(2))
scatter(lsri(1,:),lsri(2,:),mkrsz,'k','MarkerEdgeAlpha',mktrns)
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
cleanTicks
yticks([])
xlabel('{\itl}_{MB,10}');
text(0.05,0.95,'B','Units','normalized')

% Hypothetical ideal
hypI = repmat([lsri(1,:,1);lsri(2,:,1)],1,1,130);

rng(1)
noise = normrnd(0,0.001,size(hypI));
noise(2,:,:) = noise(2,:,:)/4;

hypI = hypI+noise;

axes(sp(3))
scatter(hypI(1,:),hypI(2,:),mkrsz,'k','MarkerEdgeAlpha',mktrns)
xlim([min_l_scale max_l_scale])
ylim([0 max_s_scale])
cleanTicks
yticks([])
xlabel('{\itl}_{MB,10}');
text(0.05,0.95,'C','Units','normalized')

if print_figures
    save2pdf([base,'\hypI.pdf'])
end

%%
close all

melcomp_9;

figure(1)
set(gcf, 'Position',  [plot_where plot_size.*[1,1.5]])

if print_figures
    save2pdf([base,'\ScoreByPerc.pdf'])
end

figure(2)
set(gcf, 'Position',  [plot_where plot_size.*[1,2.5]])

if print_figures
    print([base,'\OutputByPerc'],'-djpeg')
end
if print_figures
    save2pdf([base,'\OutputByPerc.pdf'])
end

%%
close all
peak_locations = melcomp_9_caller(-88:162);

figure(1)
set(gcf, 'Position',  [plot_where plot_size])

if print_figures
    save2pdf([base,'\optimality.pdf'])
end

%%

figure('Position',[plot_where plot_size.*[1,1.5]]),

subplot(1,2,1)
melcomp_9('mel_offset',516-488,'plt','single');
axis equal
axis([-2.5 2.5 -2.5 2.5])
cleanTicks

subplot(1,2,2)
melcomp_9('mel_offset',576-488,'plt','single');
axis equal
axis([-2.5 2.5 -2.5 2.5])
cleanTicks

if print_figures
    save2pdf([base,'\optimal.pdf'])
end










