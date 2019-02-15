clear, clc, close all

for i=1:12%[6:9,12]
    melcomp_2('Z_ax',i,'plt','3D');
end

%%

melcomp_2('SPD','D-series',...
    'SRF','Vrhel_nat_1',...
    'SSF','SP',...
    'plt','3D');
view(201,34)

% expectedSPD = {'Granada_sub','Granada','D-series'};
% expectedSRF = {'Vrhel_nat_1','Vrhel_nat_2','Vrhel_full','Foster'};
% expectedSSF = {'SS10','SP'};
% expectedlum = {'CIE_10','SP'};

%%
clc, clear, close all

%[pc, LMSRI] = melcomp(PF_SPD,PF_refs,PF_obs,Z_ax,plt,offset)

for i = -10:10
    offset = i*10;
    [pc(i+11), LMSRI(:,:,:,i+11)] = melcomp_2(...
        'SPD','Granada',...
        'SRF','Vrhel_nat_1',...
        'SSF','SP',...
        'Z_ax',9,...
        'mel_offset',offset);
end

%%

ex = [pc.explained];

figure,
plot(ex(3,:))