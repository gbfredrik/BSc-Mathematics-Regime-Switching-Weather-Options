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

DAT = retime(Sets(1, 1).CleanSet(:,1), 'daily', 'mean'); % TODO: Change to use DailyAverage(.)
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
temp = transpose(DAT.Degrees(21000:25839,:));

[X, FVAL] = fmincon(@(x) sum((temp - SeasonFunc(x(1), x(2), x(3), x(4), 0:4839)).^2), [18;0.005;5;0]);

%%
close ALL

figure(1)
hold on
plot(temp, 'b')
plot(SeasonFunc(X(1), X(2), X(3), X(4), 0:4839), 'r', 'LineWidth', 2)
plot(SeasonFunc(X(1), X(2), 0, 0, 0:4839))
plot(18 .* ones(4840, 1), 'g', 'LineWidth', 2)
hold off

figure(2)
hold on
subplot(2,1,1)
plot(temp - SeasonFunc(X(1), X(2), X(3), X(4), 0:4839))

subplot(2,1,2)
qqplot(temp - SeasonFunc(X(1), X(2), X(3), X(4), 0:4839))
hold off



