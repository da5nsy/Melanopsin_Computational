

figure(1)
cla reset
t_i = 1;
t_j = 1;

t_i_range = 0.2:0.005:0.6;

for i = t_i_range
    cla reset
    scatter(lsri(1,:)+i*lsri(4,:),lsri(2,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    t_i(end+1) = std(lsri(1,:)+i*lsri(4,:));
    title(i)
    drawnow
end

[t_ia,t_ib] = min(t_i);
t_i_range(t_ib)
%0.4350

t_j_range = -0.3:0.005:-0.1;

for j = t_j_range
    cla reset
    scatter(lsri(1,:),lsri(2,:)+j*lsri(4,:),[],pltc_alt(:,:)','filled','MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.6)
    t_j(end+1) = std(lsri(2,:)+j*lsri(4,:));
    title(i)
    drawnow
end

figure, plot(t_j_range,t_j(2:end))

[t_ja,t_jb] = min(t_j);
t_j_range(t_jb)
%-0.2300