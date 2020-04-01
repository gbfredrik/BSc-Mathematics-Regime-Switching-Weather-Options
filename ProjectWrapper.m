%% ProjectWrapper
%% Resetter
% Clears workspace of variables and figures on start
beep off
clearvars
close ALL

%% Settings - Run on launch
% Constants are defined through a class for easy access
settings = Settings;

%% File loading
% This section does... 
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
    Sets(1, k).FileName = sourceFiles(k);
    
    Sets(1, k).DataSet = readtable(sourceDir + sourceFiles(k));
    Sets(1, k).DataSet.Properties.VariableNames = {'Date', 'Time', 'Degrees', 'Quality'};
    Sets(1, k).DataSet.Date.Format = 'default';
    Sets(1, k).DataSet.Date = Sets(1,k).DataSet.Date + Sets(1,k).DataSet.Time;
    Sets(1, k).CleanSet = timetable(Sets(1,k).DataSet.Date, Sets(1,k).DataSet.Degrees, Sets(1,k).DataSet.Quality);
    Sets(1, k).CleanSet.Properties.VariableNames = {'Degrees', 'Quality'};
    
    Sets(1, k).ShortName = extractBefore(Sets(1, k).FileName, ' smhi-');
end

clear sourceDir sourceFiles k; % Remove variables not to be used again
%% Data parsing
% Transform the data and calculate daily average temperatures


% Set time range for source data
%Sets(1,:).TimeRange = timerange('2010-01-01', '2020-01-01');
Sets(1,1).DateStart = datetime(2010,1,1);
Sets(1,1).DateEnd = datetime(2020,1,1);
for k = 1 : length(Sets) % Iterate to parse all chosen data sets
    % TODO: Filter by dates here

    Sets(1, k).CleanSet = DailyAverage(Sets(1, k).CleanSet(:,1), settings.avgType);
end

%% Deseasoning
% Remove the seasonal component of the temperature
seasonFunction = @(a0, a1, a2, a3, t) a0 + a1 * t + a2 * sin(2 * pi / 365 * (t - a3));% + a1 * (sin(2 * pi / 365 * t)).^0.1;

X = zeros(4, length(Sets));
FVAL = zeros(1, length(Sets));

for k = 1 : length(Sets) % Iterate to deseason all chosen data sets
    [X(:, k), FVAL(1, k)] = Deseason(transpose(Sets(1, k).CleanSet.Degrees(end-3650:end,:)), seasonFunction, settings.fminconOptions);
    % TODO: Add try/catch in Deseason
end

for k = 1 : length(Sets)
    Sets(1, k).Deseasoned = Sets(1, k).CleanSet(end-3650:end,:);
    Sets(1, k).Deseasoned.Degrees = Sets(1, k).Deseasoned.Degrees - transpose(seasonFunction(X(1,k), X(2,k), X(3,k), X(4,k), 0:3650));
end

clear k
%%
close ALL

% Settings for figures
showFigures = false;
saveFigures = true;
showTref = true;
showLinTrend = true;

status = zeros(1, length(Sets));

for k = 1 : length(Sets) % Iterate to generate DAT figures
    [status(k)] = GenerateDATPlot(Sets(1, k), seasonFunction, X(:, k), ...
        showFigures, saveFigures, showTref, showLinTrend);
    %fprintf(sprintf('DAT plot status: %d.\n', status(1, k)))
end


figure()
hold on
for k = 1 : length(Sets) % Iterate to generate deseasoned figures
subplot(3,2,2*k-1)
plot(Sets(1,k).Deseasoned.Time, Sets(1,k).Deseasoned.Degrees)

subplot(3,2,2*k)
qqplot(Sets(1,k).Deseasoned.Degrees)

end
hold off

clear k showFigures saveFigures showTref showLinTrend status
%%






