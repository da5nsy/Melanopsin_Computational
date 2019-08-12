clear, clc, close all

wholeset = 1;
[~,~,~,~,MB_star1] = transformToIllIndSpace(0,wholeset);
[~,~,~,~,MB_star2] = transformToIllIndSpace(590-488,wholeset);

close all

%%

KMM = KMeansMark(MB_star1)

KMM = KMeansMark(MB_star2)