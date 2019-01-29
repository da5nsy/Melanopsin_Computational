clear, close, clc

d = dir('opt*.mat');

figure, hold on
for i = 1:length(d)
    clear out
    load(d(i).name)
    plot(out(1,:),out(2,:),'-','DisplayName',d(i).name)  
end

legend

