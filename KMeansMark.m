function [KMM,km_idx] = KMeansMark(data,k,map)

if ~exist('k','var')
    k = size(data,2);
end

if ~exist('map','var')
    map = repmat(1:size(data,2),size(data,3),1)';
end

data = data - mean(data(:,:),2);  % zero mean
data = data./std(data(:,:),[],2); % normalise by SD

km_idx = kmeans(data([1,2],:)',k,'Replicates',20);
km_r = reshape(km_idx,[size(data,2),size(data,3)]); %reshape

for i=1:k
    lu(i) =  mode(km_r(map == i)); %lookup
    pcC(i) = mean(lu(i) == km_r(map == i)); %percent correct 
end

KMM = mean(pcC);

%km_m = repmat(mode(km_r')',1,size(data,3)); %mode
%d = km_r == km_m;
%KMM = mean(d(:));

end