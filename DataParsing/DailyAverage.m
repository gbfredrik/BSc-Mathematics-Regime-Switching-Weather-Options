function [DAT] = DailyAverage(tempSeries, type)
%DAILYAVERAGE Summary of this function goes here
%   Detailed explanation goes here

if isequal(type, 'MinMax')
    tempMin = retime(tempSeries, 'daily', 'min');
    tempMax = retime(tempSeries, 'daily', 'max');
    DAT = (tempMin.Degrees + tempMax.Degrees) / 2; % Todo change to return timetable still
elseif isequal(type, 'Mean')
    DAT = retime(tempSeries, 'daily', 'mean');
end

end

