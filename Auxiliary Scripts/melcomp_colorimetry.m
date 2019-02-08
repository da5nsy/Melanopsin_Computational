function [LMSRI, lsri] = melcomp_colorimetry(T_SPD, T_SRF, T_SSF, T_lum, S_sh)

T_rad = zeros([S_sh(3),size(T_SRF,2),size(T_SPD,2)]);
LMSRI = zeros([size(T_SSF,2),size(T_SRF,2),size(T_SPD,2)]);
lsri  = zeros([4,size(T_SRF,2),size(T_SPD,2)]);
t_r   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp
t_i   = zeros([2,size(T_SRF,2),size(T_SPD,2)]); %t for temp

for i=1:size(T_SPD,2)
    T_rad(:,:,i)  = T_SRF.*T_SPD(:,i);
    LMSRI(:,:,i)  = T_SSF'*T_rad(:,:,i);
    lsri(1:2,:,i) = LMSToMacBoyn(LMSRI(1:3,:,i),T_SSF(:,1:3)',T_lum');
    t_r(:,:,i)    = LMSToMacBoyn(LMSRI([1,2,4],:,i),[T_SSF(:,1:2)';T_SSF(:,4)'],T_lum');
    t_i(:,:,i)    = LMSToMacBoyn(LMSRI([1,2,5],:,i),[T_SSF(:,1:2)';T_SSF(:,5)'],T_lum');
end
lsri(3,:,:) = t_r(2,:,:); 
lsri(4,:,:) = t_i(2,:,:); 

end