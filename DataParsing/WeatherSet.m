classdef WeatherSet
    %WEATHERSET Class for temperature sets
    %   Detailed explanation goes here
    
    properties
        FileName % Name of .csv file
        ShortName
        DateStart % Start date of data
        DateEnd % Last date of data
        DataSet % Data array
        Clean
        Deseasoned
    end
end
