function exclRef = excludeSimilarRefs(lsri_avIll,minD)

%%
% figure, hold on
% scatter(lsri_avIll(1,:),lsri_avIll(2,:))
% axis equal

%%
if ~exist('minD','var')
    minD = 0.002;
end

D = pdist(lsri_avIll(1:2,:)');
exclRef = false(size(lsri_avIll,2),1);

Ds = squareform(D);
Ds(Ds == 0) = 1;

while min(D) < minD
    Ds_temp = Ds;
    Ds_temp(exclRef,:) = 1;
    Ds_temp(:,exclRef) = 1;
    [~,refToBeExcluded] = min(min(Ds_temp));
    exclRef(refToBeExcluded) = 1;
    
    D = pdist(lsri_avIll(1:2,~exclRef)');
    
%     scatter(lsri_avIll(1,~exclRef),lsri_avIll(2,~exclRef),'r.')
%     drawnow
%     pause(1)
end

