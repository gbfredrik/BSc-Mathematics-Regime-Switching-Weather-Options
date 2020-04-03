classdef WeatherSet
    %WEATHERSET Class for temperature sets
    %   Detailed explanation goes here
    
    properties
        FileName % Name of .csv file
        ShortName
        DataSet % Data array
        DateStart % Start date of data
        DateEnd % Last date of data
        Clean
        Deseasoned
    end
end
