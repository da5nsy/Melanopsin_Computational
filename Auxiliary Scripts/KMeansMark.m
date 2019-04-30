function KMM = KMeansMark(lsri)

rng(1)
km_idx = kmeans(lsri([1,2],:)',size(lsri,2),'Replicates',50);
%pltc_alt2 = pltc_alt(:,:,1);
%figure
%scatter(lsri_gw(1,:),lsri_gw(2,:),...
%    [],pltc_alt2(:,km_idx)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
%axis equal
%title('Grey world correction')
km_r = reshape(km_idx,[size(lsri,2),size(lsri,3)]); %reshape
km_m = repmat(mode(km_r')',1,size(lsri,3)); %mode
d = km_r == km_m;
KMM = mean(d(:));

end