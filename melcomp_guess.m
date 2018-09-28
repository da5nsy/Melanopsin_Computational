function melcomp_guess(offset)

% Function for trying to guess the illuminant from the S:I ratio

try
    nargin;
catch
    clear, clc, close all
    offset = 0;
end

tic
[~, LMSRI] = melcomp(1,3,1,9,NaN,offset); % [pc, LMSRI] = melcomp(PF_SPD,PF_refs,PF_obs,Z_ax,plt,offset)
toc

%% Testing S:I ratio.

cols = jet(size(LMSRI,3));

for i=1:3:size(LMSRI,3)
    scatter(log10(LMSRI(3,1:8000:end,i)),log10(LMSRI(5,1:8000:end,i)),'MarkerEdgeColor',cols(i,:))
    %scatter(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),'filled')
    hold on
    
    fit(i,:) = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
    yfit = polyval(fit(i,:),log10(LMSRI(3,:,i)));
    legend('AutoUpdate','Off')
    plot(log10(LMSRI(3,:,i)),yfit,'Color',cols(i,:))
    legend('AutoUpdate','On')
end
legend('Location','best')
axis equal
xlabel('log(S)')
ylabel('log(I)')

% figure, hold on
% for i=1:5:size(LMSRI,3)
%     scatter(lsri(2,:,i),lsri(4,:,i),'filled')
%     hold on
%     title('s,i')
% end

figure, hold on
plot(fit(1:3:end,1),'o')
plot(fit(1:3:end,2),'o')

return

%%

cols = jet(size(LMSRI,3));

for i=round(size(LMSRI,3)/2)
    scatter(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),'.','MarkerEdgeColor',cols(i,:))
    %scatter(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),'filled')
    hold on
    
    fit = polyfit(log10(LMSRI(3,:,i)),log10(LMSRI(5,:,i)),1);
    yfit = polyval(fit,log10(LMSRI(3,:,i)));
    legend('AutoUpdate','Off')
    plot(log10(LMSRI(3,:,i)),yfit,'Color',cols(i,:))
    legend('AutoUpdate','On')
end
legend('Location','best')
axis equal
xlabel('log(S)')
ylabel('log(I)')

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

end