function [f] = GenerateStatePlot(T_d, Pr1, Pr2, threshold)
%GENERATESTATEPLOT Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 3)
   threshold = 0.5; 
end

times = 1:length(T_d);

%%
base = (Pr1 >= 1-threshold) .* T_d;
shifted = (Pr2 > threshold) .* T_d;
base(base == 0) = nan;
shifted(shifted == 0) = nan;

f = figure();
subplot(2, 1, 1);
plot(times, base, 'b.', times, shifted, 'rx')
axis([1 times(end) min(T_d)-1 max(T_d)+1]);
xlabel('Days');
ylabel('Deseasoned temperature');

sub2 = subplot(2, 1, 2);
hold on
bar(times, Pr2, 'b')
plot(times, threshold .* ones(1, length(T_d)), 'r')
axis([1 times(end) 0 1]);
xlabel('Days');
ylabel('Probability of shifted state');
hold off
%sub2.Position = sub2.Position + [0 0 0 0.5];
fprintf("Optima has %i base points and %i shifted points.\n", ...
    sum(Pr1 >= 1-threshold), sum(Pr2 >= threshold));

%%
% T_d_cum = cumsum(T_d);
% base_cum = (Pr1 >= 1-threshold) .* T_d_cum;
% shifted_cum = (Pr2 > threshold) .* T_d_cum;
% base_cum(base_cum == 0) = nan;
% shifted_cum(shifted_cum == 0) = nan;
% 
% g = figure();
% plot(times, base_cum, 'b.', times, shifted_cum, 'rx')

%%
tuple = [T_d(2:end) - T_d(1:end-1), Pr1(2:end)];
tuple = sortrows(tuple);
base_changes = (tuple(:,2) >= 1-threshold) .* (tuple(:,1));
shifted_changes = (1 - tuple(:,2) >= threshold) .* (tuple(:,1));
base_changes(base_changes == 0) = nan;
shifted_changes(shifted_changes == 0) = nan;

h = figure();
hold on
bar(times(2:end), base_changes)
bar(times(2:end), shifted_changes)
hold off
%%
pause(0.001)

end

