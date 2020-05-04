%% Generate Wiener process graph
%BPATH2  Brownian path simulation: vectorized
%Credit for code goes to: https://sites.me.ucsb.edu/~moehlis/APC591/tutorials/tutorial7/node2.html

rng('default')              % set the state of randn
T = 2; N = 365*2; dt = T/N;

dW = sqrt(dt)*randn(1,N);   % increments
W = cumsum(dW);             % cumulative sum

plot([0:dt:T],[0,W],'r-')   % plot W against t
xlabel('t','FontSize',16)
ylabel('W(t)','FontSize',16,'Rotation',0)