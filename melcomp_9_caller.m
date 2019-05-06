clear, clc, close all


mel_offset_range = [-100:5:200];%-80:20:120;
for i = 1:length(mel_offset_range)
    results(:,:,i) = melcomp_9('mel_offset',mel_offset_range(i),...
        'pcSurfRange',100);
end

%%

r = squeeze(results);

%%
figure,
plot(488+   mel_offset_range,r)