% figures for paper

%% Fig: Monotonicity_concept

x  = 0:0.01:1;
y1 = (x.^2)/max(x.^2);
y2 = (sin(x*7)+1)/2;

figure('Position',[672   395   900   400])

subplot(1,2,1), axis square, hold on

plot([0.6,0.6],[0,y1(61)],'r')
plot([0,0.6],[y1(61),y1(61)],'r')

plot(x,y1,'k')

xlabel('Input'),ylabel('Determined Output')
xticks(0:0.2:1); yticks(0:0.2:1);


subplot(1,2,2), axis square, hold on

plot([0.6,0.6],[0,0.93],'r')
plot([0,0.6],[0.033,0.033],'r:')
plot([0,0.6],[0.42,0.42],'r:')
plot([0,0.6],[0.93,0.93],'r:')

plot(y2,x,'k')
xlabel('Input'),ylabel('Undetermined Output')
xticks(0:0.2:1); yticks(0:0.2:1);

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\Monotonicity_concept','-depsc')

%% Fig: True3D

rng(4) %Generate same 'random' numbers each time
x=rand(10,1);
y=rand(10,1);

figure('Position',[672   395   900   400])

subplot(1,2,1), axis square
scatter(x,y,'k','filled');
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('x'); ylabel('y');

subplot(1,2,2), axis square
scatter(x,y./y,'k','filled');
xlim([0 1]); ylim([0 1]);
xticks(0:0.2:1); yticks(0:0.2:1);
xlabel('x'); ylabel('y/y');

%print('C:\Users\cege-user\Dropbox\UCL\Ongoing Work\Melanopsin Computational\Writing\True3D','-depsc')

