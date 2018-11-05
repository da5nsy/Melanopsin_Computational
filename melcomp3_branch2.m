%% Shortcut to this section

clear, clc, close all

load melcomp_3_fullWorkspace.mat %load everything generated melcomp_3
%load melcomp_3_correlation_results.mat %load the results of the above code run (until corr_return) for all signals
clear p_1 p_2


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
                f(i,j,ill,:) = orthogonalRegress(log10(LMSRI(i,:,ill)),log10(LMSRI(j,:,ill)));
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

plt_SymmetricalLinesDemo = 0;

if plt_SymmetricalLinesDemo
    
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
            
            xlabel(['x = ',plt_lbls{i},', y = ',plt_lbls{j}])
            ylabel('PC2')
            zlabel('PC1')
            title('m scores')
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
            
            xlabel(['x = ',num2str(i),', y = ',num2str(i)])
            ylabel('PC2')
            zlabel('PC1')
            title('c scores')
        end
    end
end


%% Plot showing segmentation by s3

plt_seg = 1;
plt_ass = 1;

MorC = 1; %1=m, 2=c

NoD  = 40; %number of divisions
cols = lines(NoD); %cols = jet(n);

block = zeros(NoD,2);
fit_t = zeros(NoD,2,5,5);

for s1i = 1:5
    for s2i = 1:5
        if s1i == s2i
            continue
        end
        sc_t = [squeeze(f(s1i,s2i,:,MorC)) pc_p.score(:,2)]; %scatter temp
        
        s3 = pc_p.score(:,1)-min(pc_p.score(:,1)); %Using PC1 as sv
        mI = max(s3);
        
        
        sc    = NaN([size(sc_t) NoD]); %scatter
        
        if plt_seg
            figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
            hold on
            set(gca, 'FontSize', 16)
            
            xlabel([plt_lbls{s1i},' against ',plt_lbls{s2i}]) %!!!!!!!!!!!!!
            ylabel('PC2')
            zlabel('Normalised PC1')
            
            axis tight
            cla
        end
        
        for i=1:NoD
            block(i,[1 2]) = [(i-1)*(mI/NoD), i*(mI/NoD)]; %compute lower and upper bounds for block
            bmi = and(s3>=(block(i,1)),s3<block(i,2)); %block membership index
            sc(bmi,:,i) = sc_t(bmi,:); %sorts values into order based on block membership
            fit_t(i,:,s1i,s2i) = orthogonalRegress(sc(~isnan(sc(:,1,i)),1,i),sc(~isnan(sc(:,2,i)),2,i)); %fit temp
            fit_t(isnan(fit_t)) = 0;
            
            if any(fit_t(i,:,s1i,s2i)) && plt_seg
                x_seg = linspace(min(sc(:,1,i)),max(sc(:,1,i)),20);
                y_seg = (fit_t(i,1,s1i,s2i) * x_seg) + fit_t(i,2,s1i,s2i);
                scatter3(sc(:,1,i),sc(:,2,i),s3,[],cols(i,:),'.');
                plot3(x_seg,y_seg,repmat(mean(block(i,:)),size(x_seg,2),1),'Color',cols(i,:));
                drawnow
            end
        end
        
        if plt_ass % Assess trends in fitted lines through data
            figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
            hold on
            set(gca, 'FontSize', 16)
            
            fit_t_idx = and(fit_t(:,1,s1i,s2i),fit_t(:,2,s1i,s2i));
            x_ass = mean(block(fit_t_idx,:),2);
            y_1 = fit_t(fit_t_idx,1,s1i,s2i);
            y_2 = fit_t(fit_t_idx,2,s1i,s2i);
            scatter(x_ass,y_1,'ro')
            scatter(x_ass,y_2,'bo')
            
            xlabel([plt_lbls{s1i},' against ',plt_lbls{s2i}])
            ylabel('Value in line of best fit')
            
            p_1([1,2],s1i,s2i) = orthogonalRegress(x_ass,y_1);
            p_2([1,2],s1i,s2i) = orthogonalRegress(x_ass,y_2);
            y_1p = polyval(p_1([1,2],s1i,s2i), x_ass);
            y_2p = polyval(p_2([1,2],s1i,s2i), x_ass);
            
            plot(x_ass,y_1p,'r:')
            plot(x_ass,y_2p,'b:')
            
            
            legend('m','c',[num2str(p_1(1,s1i,s2i)),'*x + ',num2str(p_1(2,s1i,s2i))],...
                [num2str(p_2(1,s1i,s2i)),'*x + ',num2str(p_2(2,s1i,s2i))],...
                'Location','Southwest')
            
        end
    end
end


%%

for s1i = 1:5
    for s2i = 1:5
        if s1i == s2i
            continue
        end
        
        figure('Position',[plot_where 800 800],'defaultLineLineWidth',2)
        hold on
        set(gca, 'FontSize', 16)
        
        m = p_1(1,s1i,s2i) * s3 + p_1(2,s1i,s2i);
        c = p_2(1,s1i,s2i) * s3 + p_2(2,s1i,s2i);
        
        estimatedPC2 =  m .* (squeeze(f(s1i,s2i,:,MorC))) + c;
        
        scatter3(estimatedPC2,pc_p.score(:,2),s3,'k.')
        scatter3(estimatedPC2(s3>0.5),pc_p.score((s3>0.5),2),s3(s3>0.5),'g.')
        
        xlabel(['(',num2str(p_1(1,s1i,s2i),3),' * PC1 + ',num2str(p_1(2,s1i,s2i),3),...
            ') * f(',num2str(s1i),',',num2str(s2i),...
             ') + (',num2str(p_2(1,s1i,s2i),3),' * PC1 + ',num2str(p_2(2,s1i,s2i),3),')'])
        ylabel('PC2')
        zlabel('PC1')
    end
end

%% Guesstimates game

MorC = 1;

% OK, here's the game:
% You get the cs for n surfaces under 1 illuminant.
% Tell me what the PC2 score for that illuminant is.

g_nSur=4; % number of surfaces

%rng(4);
g_ref_ind = randi(size(cs,2),g_nSur,1); %surfaces index  
g_ill_ind = randi(size(cs,2),1,1);

g_cs = cs(1:5,g_ref_ind,g_ill_ind);

for i=1:5
    for j=1:5
        if i == j
            continue
        end
        g_fits(:,i,j) = orthogonalRegress(g_cs(i,:),g_cs(j,:));
        
        
        m = p_1(1,i,j) * (pc_p.score(g_ill_ind,1)-min(pc_p.score(:,1))) + p_1(2,i,j);
        c = p_2(1,i,j) * (pc_p.score(g_ill_ind,1)-min(pc_p.score(:,1))) + p_2(2,i,j);
        
        g_estimatedPC2(i,j) =  m .* g_fits(MorC,i,j) + c;
    end
end

correct_answer = pc_p.score(g_ill_ind,2)

%%
i=1;
j=3;

figure, hold on
scatter(g_cs(i,:),g_cs(j,:))

x = linspace(min(g_cs(i,:)),max(g_cs(i,:)));
y = g_fits(1,i,j)*x+g_fits(2,i,j);
plot(x,y)

% Now I should be able to look up, for example, that m value, at this value
% of pc1, and get a pc2 estimate.


%%
correct_answer = pc_p.score(g_ill_ind,2)




