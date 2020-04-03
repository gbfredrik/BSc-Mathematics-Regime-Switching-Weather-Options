function [ld] = IsLeapDay(d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ld = false(1, length(d));

for i = 1:length(d)
    if (month(d(i)) == 2) && (day(d(i)) == 29)
        ld(i) = true;
    else
        ld(i) = false;
    end
end
end

