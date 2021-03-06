function melcomp_guess(offset)

% Function for trying to guess the illuminant from the S:I ratio

% Needs updating to call melcomp_2 correctly

try
    nargin;
catch
    clear, clc, close all
    offset = 0;
end

tic
[~, LMSRI] = melcomp_2(2,2,1,9,NaN,offset); % [pc, LMSRI] = melcomp_2(PF_SPD,PF_refs,PF_obs,Z_ax,plt,offset)
toc

% PF_SPD
% 1 = CIE D series
% 2 = Hern�ndez-Andr�s+

% PF_refs
% 1 = Vhrel+ (natural only)
% 2 = Vhrel+ (all)
% 3 = Foster+

% PF_obs
% 1 = PTB Smith-Pokorny

%% Testing S:I ratio.

cols = jet(size(LMSRI,3));
ds_r = 10; %downsample reflectances. %8000 for Foster data
ds_i = 1; %downsample illums. %3 for CIE D

for i=1:ds_i:size(LMSRI,3)
    %scatter(log10(LMSRI(3,1:ds_r:end,i)),log10(LMSRI(5,1:ds_r:end,i)),'MarkerEdgeColor',cols(i,:))
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

load melcomp_3_pc_p.mat

figure, hold on
plot3(1:length(1:ds_i:size(LMSRI,3)),fit(1:ds_i:end,1),pc_p.score(1:ds_i:size(LMSRI,3),2),'o')
plot3(1:length(1:ds_i:size(LMSRI,3)),fit(1:ds_i:end,2),pc_p.score(1:ds_i:size(LMSRI,3),2),'o')
 
% figure, hold on
% plot3(1:i,fit(1:i,1),pc_p.score(1:i,2),'.')
% plot3(1:i,fit(1:i,2),pc_p.score(1:i,2),'.')

return %if this is being called from somewhere else, end here

%% Does averaging over spatially juxtaposed points provide improvement?

figure, hold on

subplot(1,2,1)
for i=1:ds_i:size(LMSRI,3)
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