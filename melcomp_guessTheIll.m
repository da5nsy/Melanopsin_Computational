clear, clc, close all

% [pc, LMSRI] = melcomp(PF_SPD,PF_refs,PF_obs,Z_ax,plt,offset)
[pc, LMSRI] = melcomp(1,3,1,9,NaN,0);

%% Testing S:I ratio.

cols = jet(size(LMSRI,3));

figure, hold on
for i=1:size(LMSRI,3)
    scatter(log2(LMSRI(3,:,i)),log2(LMSRI(5,:,i)),'MarkerEdgeColor',cols(i,:))
    %scatter(log2(LMSRI(3,:,i)),log2(LMSRI(5,:,i)),'filled')
    hold on
    title('log(S),log(I)')
    
    fit = polyfit(log2(LMSRI(3,:,i)),log2(LMSRI(5,:,i)),1);
    yfit = polyval(fit,log2(LMSRI(3,:,i)));
    plot(log2(LMSRI(3,:,i)),yfit,'Color',cols(i,:))
end
legend('Location','best')

% figure, hold on
% for i=1:5:size(LMSRI,3)
%     scatter(lsri(2,:,i),lsri(4,:,i),'filled')
%     hold on
%     title('s,i')
% end


%% Does averaging over spatially juxtaposed points provide improvement?

figure, hold on

subplot(1,2,1)
for i=1:5:size(LMSRI,3)
    scatter(log2(LMSRI(3,:,i)),log2(LMSRI(5,:,i)),'filled')
    hold on
    title('log(S),log(I)')
end
subplot(1,2,2)
for j=1:50
    cla
    
    windowSize = j;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    LMSRI_smoothed = LMSRI;
    for i=1:size(LMSRI,1) %aka 5
        for j=1:size(LMSRI,3)
            LMSRI_smoothed(i,:,j) = filter(b,a,LMSRI(i,:,j));
        end
    end
    
    for i=1:5:size(LMSRI,3)
        scatter(log2(LMSRI_smoothed(3,:,i)),log2(LMSRI_smoothed(5,:,i)),'filled')
        hold on
        title('log(S),log(I)')
    end
    pause(0.2)
end