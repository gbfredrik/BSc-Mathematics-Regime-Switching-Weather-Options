classdef WeatherSet
    %WEATHERSET Class for temperature sets
    %   Detailed explanation goes here
    
    properties
        FileName % Name of .csv file
        ShortName
        InSample % Time range of in-sample data
        OutOfSample % Time range for out-of-sample
        DataSet % Data array
        Clean % Set of cleaned DATs
        Deseasoned % Set of deseasoned DATs
        Season_Theta % Optimal deseasoning parameters
        Season_FVal % Optimal function value for deseasoning
        
        Skewness
        Kurtosis
        
        ML_Theta
        ML_FVal
        KS
        AD
        
        EM
    end
end
