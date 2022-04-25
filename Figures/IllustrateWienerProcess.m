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

%% 4 bull/bear
clear
rng(2)
bull1 = normrnd(0.01, 0.125, 335, 1);
bear1 = normrnd(0.02, 0.35, 30, 1);
bull2 = normrnd(0.01, 0.125, 335, 1);
bear2 = normrnd(0.02, 0.35, 30, 1);
bull3 = normrnd(0.01, 0.125, 335, 1);
bear3 = normrnd(0.02, 0.35, 30, 1);
bull4 = normrnd(0.01, 0.125, 335, 1);
bear4 = normrnd(0.02, 0.35, 30, 1);
returns = [bull1; bear1; bull2; bear2; bull3; bear3; bull4; bear4];
plot(returns)
%% 2 bull/1 bear
clear
rng(3)
bull1 = normrnd(1, 0.125, 100, 1);
bear1 = normrnd(2, 0.45, 15, 1);
bull2 = normrnd(1, 0.125, 65, 1);
returns = [bull1; bear1; bull2];
plot(returns)