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

plt_relBetweenSignals = 0;

if plt_relBetweenSignals
    ill = 500;
    md = max(max(cs(1:5,:,ill))); %max data
    
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
end

%% Average correlation, across all illums

plt_c_av = 0;

if plt_c_av
    figure, 
    imagesc(c_av,[min(c_av(:)) 1])
    colormap('gray')
    colorbar('YTick',[min(c_av(:)) 1])
    axis image
    xticks(1:5); xticklabels(plt_lbls(1:5));
    yticks(1:5); yticklabels(plt_lbls(1:5));
    
    title('Average correlation, across all illums')
    
    % I don't like how the diag values are 1. Distracts from meaningful values
end

%% Calculate lines of best fit under each illum

try
    load 'melcomp_3_correlationBetweenSignalsForEachIllum.mat'
catch %takes roughly 30 seconds
    f=zeros(5,5,size(cs,3),2);
    for i=1:5
        for j = 1:5
            for ill = 1:size(cs,3)
                f(i,j,ill,:) = polyfit(log10(LMSRI(i,:,ill)),log10(LMSRI(j,:,ill)),1);
            end
        end
        disp(i)
    end
    save('melcomp_3_correlationBetweenSignalsForEachIllum.mat','f')
end

plt_firstIllumMScores = 0; % Visualise the first illum m scores
if plt_firstIllumMScores 
    figure, 
    imagesc(f(:,:,1,1))
    colorbar
    axis image
    colormap('gray')
    
    xticks(1:5); xticklabels(plt_lbls(1:5));
    yticks(1:5); yticklabels(plt_lbls(1:5));
end

%%

plt_nonSymmetricalLinesDemo = 0;

if plt_nonSymmetricalLinesDemo
    
    ill = 500;
    cols = lines(2);
    
    figure, hold on
    scatter(log10(LMSRI(3,:,ill)),log10(LMSRI(2,:,ill)),...
        [],cols(1,:),'filled','MarkerFaceAlpha',0.5)
    x = linspace(min(log10(LMSRI(3,:,ill))),max(log10(LMSRI(3,:,ill))));
    y = x * f(3,2,ill,1) + f(3,2,ill,2);
    plot(x,y,'Color',cols(1,:))
    
    scatter(log10(LMSRI(2,:,ill)),log10(LMSRI(3,:,ill)),...
        [],cols(2,:),'filled','MarkerFaceAlpha',0.5)
    x = linspace(min(log10(LMSRI(2,:,ill))),max(log10(LMSRI(2,:,ill))));
    y = x * f(2,3,ill,1) + f(2,3,ill,2);
    plot(x,y,'Color',cols(2,:))
    
    plot(x,x,'k:')
    
    axis image
end

% %%
%
% polyfit(log10(LMSRI(2,:,ill)),log10(LMSRI(3,:,ill)),1)
% polyfit(log10(LMSRI(3,:,ill)),log10(LMSRI(2,:,ill)),1)
%
% %%
% clear, close all
%
% rng(1)
% x = linspace(0,1,25);
% y = rand(1,25);
%
% figure, hold on
% scatter(x,y)
% scatter(y,x)
%
% p1 = polyfit(x,y,1);
% p2 = polyfit(y,x,1);
%
% plot(x,x*p1(1)+p1(2))
% plot(y,y*p2(1)+p2(2))
%
% axis image
% legend
%
% plot(x,x,'k:')
%
%
% %%
%
% clear, close all
%
% x = linspace(0,1,25);
% y = x*2-1;
% y(5) = -0.9;
% y(2) = -0.5;
% y(4) = -0.2;
% y(9) = -0.7;
%
% figure, hold on
% scatter(x,y)
% scatter(y,x)
%
% p1 = polyfit(x,y,1);
% p2 = polyfit(y,x,1);
%
% plot(x,x*p1(1)+p1(2))
% plot(y,y*p2(1)+p2(2))
%
% axis image
% legend('Location','best')
%
% plot(x,x,'k:')

%% Plot the spiral for each m and c score

plt_m_scores = 0;
plt_c_scores = 0;

if plt_m_scores
    for i=1:5
        for j=1:5
            if i==j
                continue
            end
            figure,
            scatter3(f(i,j,:,1),pc_p.score(:,2),pc_p.score(:,1),'k.')
            
            xlabel(plt_lbls(i))
            ylabel(plt_lbls(j))
            zlabel('PC1')
        end
    end
end

if plt_c_scores
    for i=1:5
        for j=1:5
            if i==j
                continue
            end
            figure,
            scatter3(f(i,j,:,2),pc_p.score(:,2),pc_p.score(:,1),'k.')
            
            xlabel(plt_lbls(i))
            ylabel(plt_lbls(j))
            zlabel('PC1')
        end
    end
end
