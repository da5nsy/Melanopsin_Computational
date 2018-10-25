%%
clear
clc
close all

for i=1:21
    %for j=1:5
    j=1;
        
        close all
        %clc
        
        [corr_return(i,j), p_1(i,j,:), p_2(i,j,:)] = melcomp_3(i,j);
    %end
end

%%
plt_lbls{1}  = 'L'; %writing out this way so that there's a quick reference as to which value of Z_ax does what
plt_lbls{2}  = 'M';
plt_lbls{3}  = 'S';
plt_lbls{4}  = 'R';
plt_lbls{5}  = 'I';
plt_lbls{6}  = 'l';
plt_lbls{7}  = 's';
plt_lbls{8}  = 'r';
plt_lbls{9}  = 'i';
plt_lbls{10} = 'L+M';
plt_lbls{11} = '(0.6373*L)+(0.3924*M)';
plt_lbls{12} = 'r + i';
plt_lbls{13} = 'L/M';
plt_lbls{14} = 'L/S';
plt_lbls{15} = 'L/R';
plt_lbls{16} = 'L/I';
plt_lbls{17} = 'S/R';
plt_lbls{18} = 'S/I';
plt_lbls{19} = '(S+I)/(L+M)';
plt_lbls{20} = 'S/I fit m';
plt_lbls{21} = 'S/I fit c';

figure,
plot(abs(corr_return))
ylim([0 1])
ylabel('Correlation')
xticks(1:21)
xticklabels(plt_lbls)
set(gca,'xaxisLocation','top')
xtickangle(45)

legend(plt_lbls{1:5})

%%

save('melcomp_3_correlation_results')