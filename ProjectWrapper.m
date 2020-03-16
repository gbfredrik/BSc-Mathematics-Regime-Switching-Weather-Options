%% Resetter
% Empties workspace of variables on start
clearvars

%% File loading
% This section does... wanka
sourceDir = "DataSets/";
sourceFiles = ["Stockholm Bromma smhi-opendata 20200130.csv", ...
               "Kiruna Flygplats smhi-opendata 20200130.csv", ...
               "Falsterbo smhi-opendata 20200130.csv"];

if (~exist(sourceDir, 'dir'))
    fprintf("DataSets dir not found!\n\n")
    return
end

Sets(1:length(sourceFiles)) = WeatherSet;
for k = 1 : length(Sets) % Iterate to load all chosen data sets
    %TODO load date: Sets(1, k).DateStart = readcell(sourceFiles(k), 'range', 'A8');
    %TODO load date: Sets(1, k).DateEnd = xlsread(sourceFiles(k), 'B8:B8');
    
    Sets(1, k).DataSet = readtable(sourceDir + sourceFiles(k));
    Sets(1, k).DataSet.Properties.VariableNames = {'Date', 'Time', 'Degrees', 'Quality'};
    Sets(1, k).DataSet.Date.Format = 'default';
    Sets(1, k).DataSet.Date = Sets(1,k).DataSet.Date + Sets(1,k).DataSet.Time;
    Sets(1, k).CleanSet = timetable(Sets(1,k).DataSet.Date, Sets(1,k).DataSet.Degrees, Sets(1,k).DataSet.Quality);
    Sets(1, k).CleanSet.Properties.VariableNames = {'Degrees', 'Quality'};
    
    Sets(1, k).FileName = sourceFiles(k);
end
clear sourceDir sourceFiles k; % Remove variables not to be used again
%% Data parsing
% Transform the data and calculate daily average temperatures

% Set time range for source data
%Sets(1,:).TimeRange = timerange('2000-01-01', '2020-01-01');
Sets(1,1).DateStart = datetime(2000,1,1);
Sets(1,1).DateEnd = datetime(2020,1,1);

%% Deseasoning
% Remove the seasonal component of the temperature
SeasonFunc = @(a0, a1, a2, a3, t) a0 + a1 * t + a2 * sin(2 * pi / 365 * (t - a3));% + a1 * (sin(2 * pi / 365 * t)).^0.1;
%figure(1);
%plot(SeasonFunc(18, 0.000005, 0.5, 0.5, 1:0.1:1000))

[x, FVAL] = fmincon(@(x) (Sets(1,:).TimeRange - SeasonFunc(x(1), x(2), x(3), x(4), 1:365*20))^2, [18;0.005;5;0], [], [], [])

%%

