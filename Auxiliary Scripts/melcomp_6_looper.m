function  pc = melcomp_6_looper(offset,norm,plt) 

if ~exist('offset','var') %do this properly with nargin !!!!!!!!!!!!!
    disp('No offset passed. Setting offset to 0.')
    offset = 0;
    disp('No norm command passed. Setting norm to 1 (positive)')
    norm = 1;
end

%% Data

[T_SPD, T_SRF, T_SSF, S_sh] = melcomp_loader('SPD','Granada','SRF','Vrhel_nat_2','SSF','SS10','mel_offset',offset);
T_SPD = T_SPD(:,1:20:end);

%

sf_10 = [0.69283932, 0.34967567, 0.05547858]; %energy 10deg from CIE 170-2:2015
T_lum = sf_10(1)*T_SSF(:,1)+sf_10(2)*T_SSF(:,2);

%

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_SSF,2),size(T_SRF,2),size(T_SPD,2)]);
lsri  = zeros([4,size(T_SRF,2),size(T_SPD,2)]);
t_r   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp
t_i   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_SSF'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoynDG(LMSRI(1:3,:,i),T_SSF(:,1:3)',T_lum');
    t_r(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
    t_i(:,:,i)    = LMSToMacBoynDG(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
end
lsri(3,:,:) = t_r(2,:,:); clear t_r
lsri(4,:,:) = t_i(2,:,:); clear t_i

%% correct axes for variance

% figure,
% scatter3(lsri(1,:),lsri(2,:),lsri(4,:));

lsri = lsri(:,:);

if norm
    for i=1:4
        lsri(i,:) = (lsri(i,:) - mean(lsri(i,:)))./std(lsri(i,:));
    end
end


%% pca

lsi = lsri([1,2,4],:);

[pc.coeff, pc.score, pc.latent, pc.tsquared, pc.explained, pc.mu] = pca(lsi');

%%

if plt
    %figure,
    scatter3(lsri(1,:),lsri(2,:),lsri(4,:),'k','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6);
end
    

end