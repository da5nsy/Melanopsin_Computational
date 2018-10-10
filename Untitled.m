
melcomp_3;
close all
%%

SoI = squeeze(mean(cs(18,:,:),2));
SoIL = log2(SoI);
pc1 = pc.SCORE(:,1);
pc2 = pc.SCORE(:,2);
pc1n = (pc1-min(pc1)).^(1/4);

% figure, hold on
% scatter3(SoI,pc2,pc1,'r.')
% axis tight
% xlabel('SoI')
% ylabel('PC2')
% zlabel('PC1')
%
% figure, hold on
% scatter3(SoI,pc2,pc1n,'b.')
% axis tight
% xlabel('SoI')
% ylabel('PC2')
% zlabel('PC1n')
% axis equal

figure, hold on
scatter3(SoIL,pc2,pc1n,'b.')
axis tight
xlabel('SoIL')
ylabel('PC2')
zlabel('PC1n')
axis equal

%% Segment
close all
clc

n = 40;
m1 = max(pc1n);
cols = flag(n); %cols = jet(n);

sc_t = [SoIL pc2]; %scatter temp
sc = NaN([size(sc_t) n]);
%sc = repmat(sc_t,1,1,n);

figure, hold on
axis equal

for i=1:n
    sc(and(pc1n<(i*(m1/n)),pc1n>((i-1)*(m1/n))),:,i) = sc_t(and(pc1n<(i*(m1/n)),pc1n>((i-1)*(m1/n))),:);
    fit_t(i,:) = polyfit(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i),1); %fitty, lol. Seriously though, temp.
    
    x = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
    y = (fit_t(i,1) * x) + fit_t(i,2);
    
    if any(fit_t(i,:))
        scatter(sc(:,1,i),sc(:,2,i),[],cols(i,:));
        
        plot(x,y,'Color',cols(i,:));
    end
    drawnow
    pause(0.2)
end

%%

figure, hold on
plot(5:39,fit_t(5:39,1))
plot(5:39,fit_t(5:39,2))

a = polyfit(5:39,fit_t(5:39,1)',6);
y_p = polyval(a,5:39);
plot(5:39,y_p);




