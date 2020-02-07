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
    Sets(1, k).FileName = sourceFiles(k);
end
clear sourceDir sourceFiles k; % Remove variables not to be used again
%% Data parsing
% Transform the data and calculate daily average temperatures


%% Deseasoning
% Remove the seasonal component of the temperature
SeasonFunc = @(a0, a1, a2, a3, t) a0 * sin(2 * pi / 365 * (t - a1)) + a2*t + a3;% + a1 * (sin(2 * pi / 365 * t)).^0.1;
figure(1);
plot(SeasonFunc(14, 1, 0.00005, 8, 1:0.1:10000))

%%

