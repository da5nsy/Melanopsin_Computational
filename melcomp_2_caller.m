clear, clc, close all

for i=1:12%[6:9,12]
    melcomp_2(i)
end

%%

melcomp_2(1,1,1,1)
view(75,20)

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

%%
clc, clear, close all

%[pc, LMSRI] = melcomp(PF_SPD,PF_refs,PF_obs,Z_ax,plt,offset)

for i = -10:10
    offset = i*10;
    [pc(i+11), LMSRI(:,:,:,+11)] = melcomp_2(2,1,1,9,0,offset);
end

%%

ex = [pc.explained];

figure,
plot(ex(3,:))