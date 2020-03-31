function [DAT] = DailyAverage(tempSeries, type)
% DAILYAVERAGE Summary of this function goes here
%   Detailed explanation goes here

if isequal(type, 'MinMax')
    tempMin = retime(tempSeries, 'daily', 'min');
    tempMax = retime(tempSeries, 'daily', 'max');
    tempMin.Degrees = (tempMin.Degrees + tempMax.Degrees) / 2;
    DAT = tempMin;
elseif isequal(type, 'Mean')
    DAT = retime(tempSeries, 'daily', 'mean');
else
    DAT = 0;
    fprintf(strcat('There is no support for the type ', type, '.\n'))
    fprintf('Try "MinMax" or "Mean" instead.\n')
end
end
