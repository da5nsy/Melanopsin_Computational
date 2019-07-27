function LOCS = melcomp_9_caller(mel_offset_range)

if ~exist('mel_offset_range','var')
    mel_offset_range = [-100:200];
end

try load('melcomp_9_caller_results5')
catch
    for i = 1:length(mel_offset_range)
        [results(:,:,i),~,~,sf_s(i),sf_l(i)] = melcomp_9('mel_offset',mel_offset_range(i),...
            'pcSurfRange',100,'plt','none');
        disp(mel_offset_range(i))
    end
    r = squeeze(results);
    save('melcomp_9_caller_results5','r');
end

set(0,'defaultAxesFontName', 'Courier')

%%

hold on
plot(488+mel_offset_range,r(1,:),'Color',[0.5,0.5,0.5],'LineWidth',2) %488 max spec sens of melanopsin
plot(488+mel_offset_range,r(4,:),'k','LineWidth',2) 
axis tight

load T_cones_ss10.mat
S_x = SToWls(S_cones_ss10);
[~,max_idx] = max(T_cones_ss10');
peaks = S_x(max_idx);

for i=1:3
    plot([peaks(i),peaks(i)],[min(ylim),max(ylim)],'k:')
end
plot([488,488],[min(ylim),max(ylim)],'k--')

xlabel('Wavelength (nm)')
ylabel('k-means-mark')

% figure,
% plot(488+mel_offset_range,sf_s)
% title('sf_s')
% 
% figure,
% plot(488+mel_offset_range,sf_l)
% title('sf_l')
%% Find peaks

[~,LOCS] = findpeaks(r(4,:),488+mel_offset_range,'MinPeakHeight',0.6);

end