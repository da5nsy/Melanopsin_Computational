clear, clc, close all

wholeset = 1;
[~,~,~,~,MB_star1] = transformToIllIndSpace(0,wholeset);
[~,~,~,~,MB_star2] = transformToIllIndSpace(590-488,wholeset);

%%


scatter(MB_star1(1,:),MB_star1(2,:))