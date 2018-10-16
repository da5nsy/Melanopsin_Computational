function Colorimetric_effect_of_PC_shifting(pc_p)

% Need to make a decision on which PCA to use as input here.
% This was extracted from melcomp_OverviewPresentationFigs, where the SPD
% had been cut down to 17:97 at the first opportunity

plot_where = [20,60];

%% Show colorimetric effect of principal components

%load CMF
load T_xyz1931.mat

%start figure
figure('Position',[plot_where 800 800]), hold on
set(gca, 'FontSize', 16)
set(gcf,'defaultLineLineWidth',2)
axis square
xlim([0 1]); xticks([0,1]);
ylim([0 1]); yticks([0,1]);
xlabel('x')
ylabel('y')

%plot spectral locus
plot(T_xyz1931(1,:)./sum(T_xyz1931),T_xyz1931(2,:)./sum(T_xyz1931),...
    'ko-','MarkerFaceColor','k','MarkerSize',2);

origin.score = median(pc_p.score(:,1:3));
origin.spd   = origin.score * pc_p.coeff(:,1:3)' + pc_p.mu; %plot(origin.spd)
origin.XYZ   = T_xyz1931 * origin.spd';
origin.xy    = [origin.XYZ(1)./sum(origin.XYZ);origin.XYZ(2)./sum(origin.XYZ)];
%scatter(origin.xy(1),origin.xy(2),'k*')

for i=1:3
    n.incpc(:,i) = [linspace(min(pc_p.score(:,i)),median(pc_p.score(:,i)),20),...
        linspace(median(pc_p.score(:,i)),max(pc_p.score(:,i)),20)];
    % - incremental pc weights
    % - Should find a way to remove middle repeating value
    % - 5th/95th percentile might be preferable to min/max
end
%figure, plot(n.incpc)

n.score = repmat(origin.score,length(n.incpc),1,3);
%n.score = zeros(40,3,3);
for i=1:3
    n.score(:,i,i) = n.incpc(:,i);
end

n.spd = zeros(length(n.incpc),3,size(pc_p.score,2));

for i=1:length(n.incpc)
    for j=1:3
        n.spd(i,j,:) = n.score(i,:,j) * pc_p.coeff(:,1:3)' + pc_p.mu;
    end
end

% figure, hold on, axis tight
% plot(SToWls(S_SPD),squeeze(n.spd(:,1,:))'/max(max(squeeze(n.spd(:,1,:)))),'k')
% plot(SToWls(S_SPD),squeeze(n.spd(:,2,:))'/max(max(squeeze(n.spd(:,2,:)))),'b')
% plot(SToWls(S_SPD),squeeze(n.spd(:,3,:))'/max(max(squeeze(n.spd(:,3,:)))),'r')

for i=1:length(n.incpc)
    for j=1:3
        n.XYZ(i,j,:) = T_xyz1931 * squeeze(n.spd(i,j,:));
        n.xy(i,j,:)  = [n.XYZ(i,j,1)./sum(n.XYZ(i,j,:));n.XYZ(i,j,2)./sum(n.XYZ(i,j,:))];
    end
end

scatter3(n.xy(:,1,1),n.xy(:,1,2),n.XYZ(:,1,2),'k*')
zlabel('Y')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

scatter3(n.xy(:,2,1),n.xy(:,2,2),n.XYZ(:,2,2),'b*') %could be nice to use 'real' colours here(?)
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure

scatter3(n.xy(:,3,1),n.xy(:,3,2),n.XYZ(:,3,2),'r*')
if p, print([base,'\',num2str(p)],ff); p=p+1; end %save figure
end