figure, hold on

for i=-50:5:50
    clf
    melcomp_guess(i)
    drawnow
    disp(488 + i)
end