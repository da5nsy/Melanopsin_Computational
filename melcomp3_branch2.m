%% Shortcut to this section

clear, clc, close all

load melcomp_3_fullWorkspace.mat %load everything generated melcomp_3
%load melcomp_3_correlation_results.mat %load the results of the above code run (until corr_return) for all signals

%% Correlation matrix between the different sigansl

c_each = zeros(5,5,size(cs,3));
for i=1:size(cs,3)
    c_each(:,:,i) = corr(cs(1:5,:,i)');
end
c_av = mean(c_each,3); %mean across all illuminants

%% Plot the different signals against eachother and look at correlation (FOR ILLUMINANT #1)

ill = 500;

md = max(max(cs(1:5,:,ill))); %max data

%figure, scatter(cs(1,range,1),cs(2,range,1))

figure('Name',['Correlation between signals for illuminant #',num2str(ill)],'Color','white')
for j=1:5
    for i = 1:5
        if i==j
            continue
        end
        subplot(5,5,i+(j-1)*5)
        scatter(cs(j,:,ill),cs(i,:,ill),'k.') %only plotting points for 1 illuminant
        xlabel(plt_lbls(j))
        ylabel(plt_lbls(i))
        
        
        xlim([0 md])
        ylim([0 md])
        xticklabels([])
        yticklabels([])        
        set(gca,'Color',repmat(c_each(i,j,ill),3,1))
        
        text(md/10,md-md/10,num2str(c_each(i,j,ill),2))
    end
end

%% Average correlation, across all illums

figure, imagesc(c_av,[min(c_av(:)) 1])
colormap('gray')
colorbar('YTick',[min(c_av(:)) 1])
axis image
xticks(1:5); xticklabels(plt_lbls(1:5));
yticks(1:5); yticklabels(plt_lbls(1:5));

title('Average correlation, across all illums')

% I don't like how the diag values are 1. Distracts from meaningful values

%%

for i=1:size(LMSRI,3)
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
end 




