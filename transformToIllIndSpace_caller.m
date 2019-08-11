clc, clear, close all

%%
range= [-150:10:-60,-58:2:138,140:10:200];
wholeset = 1;

MB1_minSD = zeros(1,length(range));
MB2_minSD = zeros(1,length(range));

for i= 1:length(range)
    [MB1_minSD(i),MB2_minSD(i)]=transformToIllIndSpace(range(i),wholeset,0,0);
    disp(range(i))
end

%%
figure, hold on
plot(488+range,MB1_minSD,'.')
plot(488+range,MB2_minSD,'.')

xlabel('Peak sensitivity (nm)')
ylabel('Min SD')
yticks(ylim)
legend('k1','k2')

%%

if wholeset
    save2pdf('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\transformToIllIndSpace\offsetrange1.pdf')
else
    save2pdf('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\transformToIllIndSpace\offsetrange2.pdf')
end

%% Test out the proposed values

[~,minloc1] = min(MB1_minSD);
transformToIllIndSpace(minloc1,wholeset,1,0);
if wholeset
    save2pdf('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\transformToIllIndSpace\correctedChromaticities_range592.pdf')
end

[~,minloc2] = min(MB2_minSD);
transformToIllIndSpace(minloc2,wholeset,1,0);
if wholeset
    save2pdf('C:\Users\cege-user\Dropbox\Documents\MATLAB\Melanopsin_Computational\figs\transformToIllIndSpace\correctedChromaticities_range472.pdf')
end