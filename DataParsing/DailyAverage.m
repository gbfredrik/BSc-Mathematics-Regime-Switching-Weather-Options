function [DAT] = DailyAverage(tempSeries, type)
%DAILYAVERAGE Summary of this function goes here
%   Detailed explanation goes here

if isequal(type, 'MinMax')
   DAT = (retime(tempSeries, 'daily', 'min') + retime(tempSeries, 'daily', 'max')) / 2;
elseif isequal(type, 'Mean')
   DAT = retime(tempSeries, 'daily', 'mean');
end

end

