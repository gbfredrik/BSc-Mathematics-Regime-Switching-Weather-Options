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
sourceFiles = ["Stockholm Bromma smhi-opendata 20200404.csv", ...
               "Kiruna Flygplats smhi-opendata 20200404.csv", ...
               "Falsterbo smhi-opendata 20200404.csv"];

if (~exist(sourceDir, 'dir'))
    fprintf("DataSets dir not found!\n\n")
    return
end

Sets(1:length(sourceFiles)) = WeatherSet;
for k = 1 : length(Sets) % Iterate to load all chosen data sets
    Sets(1,k).FileName = sourceFiles(k);
    
    Sets(1,k).DataSet = readtable(sourceDir + sourceFiles(k));
    Sets(1,k).DataSet.Properties.VariableNames = {'Date', 'Time', 'Degrees', 'Quality'};
    Sets(1,k).DataSet.Date.Format = 'default';
    Sets(1,k).DataSet.Date = Sets(1,k).DataSet.Date + Sets(1,k).DataSet.Time;
    Sets(1,k).DataSet.Time = []; % Remove now unneeded time column
    
    Sets(1,k).DataSet = Sets(1,k).DataSet(Sets(1,k).DataSet.Date >= ...
        datetime(2005,01,01),:);
    
    Sets(1,k).Clean = timetable(Sets(1,k).DataSet.Date, ...
        Sets(1,k).DataSet.Degrees, Sets(1,k).DataSet.Quality);
    Sets(1,k).Clean.Properties.VariableNames = {'Degrees', 'Quality'};
    
    Sets(1,k).ShortName = extractBefore(Sets(1,k).FileName, ' smhi-');
end

clear sourceDir sourceFiles k; % Remove variables not to be used again
%% Data parsing
% Transform the data and calculate daily average temperatures

for k = 1 : length(Sets) % Iterate to parse all chosen data sets
    Sets(1,k).InSample = datetime(2008,01,01):datetime(2014,12,31);
    Sets(1,k).OutOfSample = datetime(2015,01,01):datetime(2019,12,31);
    Sets(1,k).InSample(month(Sets(1,k).InSample) == 2 ...
        & day(Sets(1,k).InSample) == 29) ...
        = []; % Clean leap days
    Sets(1,k).OutOfSample(month(Sets(1,k).OutOfSample) == 2 ...
        & day(Sets(1,k).OutOfSample) == 29) ...
        = []; % Clean leap days
    
    Sets(1,k).Clean = DailyAverage(Sets(1,k).Clean(:,1), settings.avgType);
    Sets(1,k).Clean(month(Sets(1,k).Clean.Time) == 2 ...
        & day(Sets(1,k).Clean.Time) == 29,:) ...
        = [];
end

% Check for nans
for k = 1 : length(Sets)
    n = sum(isnan(Sets(1,k).Clean.Degrees));
    if n ~= 0
        fprintf(2, sprintf('Set %s contains %d NaN values.\n', ...
            Sets(1,k).ShortName, n))
    end
end

clear k n
%% Deseasoning
% Remove the seasonal component of the temperature
seasonFunction = @(a, t) a(1) + a(2) * t + ...
    a(3) * sin(2 * pi / 365 * (t - a(4)));

guess = [8, 0.0005, 5, 100];
for k = 1 : length(Sets) % Iterate to deseason all chosen data sets
    [Sets(1,k).Season_Theta, Sets(1,k).Season_FVal] = ...
        Deseason(transpose(Sets(1,k).Clean.Degrees(Sets(1,k).InSample)), ...
        seasonFunction, length(Sets(1,k).InSample), guess, ...
        settings.fminconOptions);
end

for k = 1 : length(Sets)
    Sets(1,k).Deseasoned = Sets(1,k).Clean(Sets(1,k).InSample,:);
    Sets(1,k).Deseasoned.Degrees = Sets(1,k).Deseasoned.Degrees - ...
        transpose(seasonFunction(Sets(1,k).Season_Theta, 0:length(Sets(1,k).InSample)-1));
end

clear k guess
%% Generation of DAT and Deseasoned plots
% Allows for specific settings for the plots. 
% Used to generate figures for the written thesis report.
close ALL

% Settings for figures
showFigures = false;
saveFigures = false;
showSeason = true;
showTref = false;
showLinTrend = true;
setPeriod = "In"; % Alternatives: "In", "InOut"

status = zeros(1,length(Sets));
for k = 1 : length(Sets) % Iterate to generate DAT figures
    [status(k)] = GenerateDATPlot(Sets(1,k), seasonFunction, Sets(1,k).Season_Theta, ...
        showFigures, saveFigures, showSeason, showTref, showLinTrend, ...
        setPeriod);
    %fprintf(sprintf('DAT plot status: %d.\n', status(1,k)))
end

for k = 1 : length(Sets) % Iterate to generate deseasoned figures
    [status(k)] = GenerateDeseasonedPlots(Sets(1,k), ...
        showFigures, saveFigures);
    %fprintf(sprintf('Deseasoned plot status: %d.\n', status(1,k)))
end

for k = 1 : length(Sets) % Iterate to generate skewness and kurtosis
    Sets(1,k).Skewness = skewness(Sets(1,k).Deseasoned.Degrees);
    Sets(1,k).Kurtosis = kurtosis(Sets(1,k).Deseasoned.Degrees);
    fprintf(sprintf('Set %i ''%s'', skewness: %d, kurtosis: %d.\n', k, Sets(1,k).ShortName, Sets(1,k).Skewness, Sets(1,k).Kurtosis));
end

clear k showFigures saveFigures showSeason showTref showLinTrend ...
    setPeriod status
%% Parameter Estimation
% Define GH probability density function
f_GH = @(x, lambda, alpha, beta, delta, mu) ...
    sqrt(alpha^2 - beta^2)^(lambda) ... % TODO: SKA DET VARA SQRT HÄR?
    / (sqrt(2*pi) * delta^lambda * alpha^(lambda - 1/2) * besselk(lambda, delta * sqrt(alpha^2 - beta^2))) ...
    * (sqrt(delta^2 + (x - mu)^2))^(lambda - 1/2) ...
    * besselk(lambda - 1/2, alpha * sqrt(delta^2 + (x - mu)^2)) * exp(beta * (x - mu));

f_HYP = @(x, alpha, beta, delta, mu) f_GH(x, 1, alpha, beta, delta, mu);
f_NIG = @(x, alpha, beta, delta, mu) f_GH(x, -1/2, alpha, beta, delta, mu);

% Alternative VG pdf
f_VG = @(x, lambda, alpha, beta, mu) ...
    sqrt(alpha^2 - beta^2)^(2 * lambda) ...
    * abs(x - mu)^(lambda - 1/2) * besselk(lambda - 1/2, alpha * abs(x - mu)) ...
    / (sqrt(pi) * gamma(lambda) * (2 * alpha)^(lambda - 1/2)) ...
    * exp(beta * (x - mu));

%% Maximum likelihood estimation of 
% lambda, alpha, beta, delta, mu
guess.GH = [1, 1.7178, 0.3921, 1.6783, 0.6179]';
guess.HYP = [(1), 1.7178, 0.3921, 1.6783, 0.6179]';
guess.NIG = [(-1/2), 1.7178, 0.3921, 1.6783, 0.6179]';
guess.VG = [10, 2, -0.5, 5]';

for k = 1 : length(Sets)
     %   A, b, Aeq, beq, lb, ub, nonlcon, options
    [Sets(1, k).ML_Theta.GH, Sets(1, k).ML_FVal.GH] = fmincon(@(x) ...
        -log_likelihood_f_GH(Sets(1, k).Deseasoned.Degrees, x(1), x(2), x(3), x(4), x(5)), ...
        guess.GH, ...
        [0, -1, -1, 0 0; 0, -1, 1, 0, 0], [0;0], [], [], [-inf -inf -inf 0 -inf], [], [], settings.fminconOptions);
    Sets(1, k).ML_FVal.GH = -Sets(1, k).ML_FVal.GH; % Correct FVal to not be minimized
    
    [Sets(1, k).ML_Theta.HYP, Sets(1, k).ML_FVal.HYP] = fmincon(@(x) ...
        -log_likelihood_f_GH(Sets(1, k).Deseasoned.Degrees, x(1), x(2), x(3), x(4), x(5)), ...
        guess.HYP, ...
        [0, -1, -1, 0 0; 0, -1, 1, 0, 0], [0;0], [], [], [1 -inf -inf -inf -inf], [1 inf inf inf inf], [], settings.fminconOptions);
    Sets(1, k).ML_FVal.HYP = -Sets(1, k).ML_FVal.HYP;

    [Sets(1, k).ML_Theta.NIG, Sets(1, k).ML_FVal.NIG] = fmincon(@(x) ...
        -log_likelihood_f_GH(Sets(1, k).Deseasoned.Degrees, x(1), x(2), x(3), x(4), x(5)), ...
        guess.NIG, ...
        [0, -1, -1, 0 0; 0, -1, 1, 0, 0], [0;0], [], [], [-1/2 -inf -inf -inf -inf], [-1/2 inf inf inf inf], [], settings.fminconOptions);
    Sets(1, k).ML_FVal.NIG = -Sets(1, k).ML_FVal.NIG;
%     
    [Sets(1, k).ML_Theta.VG, Sets(1, k).ML_FVal.VG] = fmincon(@(x) ...
         -log_likelihood_f_VG(Sets(1, k).Deseasoned.Degrees, x(1), x(2), x(3), x(4)), ...
         guess.VG, ...
         [0, -1, -1, 0; 0, -1, 1, 0], [0;0], [], [], [0 0 -inf -inf], [], [], settings.fminconOptions);
     Sets(1, k).ML_FVal.VG = -Sets(1, k).ML_FVal.VG;
end
%%
%{
clearvars Pr1 Pr2 Pr1T Pr2T Theta_f p Q iter_f
model = 1;
useKim = 0;

if (model == 1)
    Theta = [0, 0.5, 0.98, 20, 2, 0.6;
             0, 0.5, 0.98, 20, 2, 0.6;
             0, 0.5, 0.98, 20, 2, 0.6]; % Initial parameter guess. Use last known 
elseif (model == 2) 
    Theta = [1, -5, 0.5, 0.98, 20, 3, 0.6;
             1, -5, 0.5, 0.98, 20, 3, 0.6;
             1, -5, 0.5, 0.98, 20, 3, 0.6];
end

p = [0.99 0.01;
     0.70 0.30];
%iter_f = zeros(1, 3);

for k = 1%length(Sets)
    [Pr1(:,k), Pr2(:,k), Pr1T(:,k), Pr2T(:,k), Theta_f, p, Q(:,k), iter_f(1,k)] = ...
        EM( ...
            Sets(1,k).Deseasoned.Degrees, ...
            Theta(k,:), ...
            p, ...
            model, ...
            100, ...
            true, ...
            useKim);
end
%}
%%
%{
% ML_Theta:
m = 1;
ml_t = Sets(1,1).ML_Theta.GH;

% GH/HYP/NIG
lambda = ml_t(1);
alpha = ml_t(2);
beta = ml_t(3);
delta = ml_t(4);
mu = ml_t(5);


% VG
% lambda = ml_t(1);
% alpha = ml_t(2);
% beta = ml_t(3);
% mu = ml_t(4);

%f_GH(1, lambda, alpha, beta, delta, mu)
pdf_series = zeros(length(Sets(1,1).Deseasoned.Degrees),1);
sorted_deg = sort(Sets(1,1).Deseasoned.Degrees);
%for i = 1:length(sorted_deg)
n = 1;
test = [];

for x_i = sorted_deg'
    %data = [ data f_GH(x_i, lambda, alpha, beta, delta, mu)];
    %pdf_series(i,1) = f_GH(sorted_deg(i), lambda, alpha, beta, delta, mu);
    %test(n,1) = f_GH(i, lambda, alpha, beta, 0, mu);
    test(n,1) = f_GH(x_i, lambda, alpha, beta, delta, mu);
    %test(n,1) = f_VG(i, lambda, alpha, beta, mu);
    n = n + 1;
end

figure()
qqplot(test)
%}
%%



%% Testing section
%{
k = 1;
min(Sets(1,k).Deseasoned.Degrees)
max(Sets(1,k).Deseasoned.Degrees)
T_deseason = Sets(1,k).Deseasoned.Degrees;

%

for t = 2 : length(T_deseason)
    kappa = Theta_f(1);
    sigma_1 = Theta_f(2);
    mu_2 = Theta_f(4);
    sigma_2 = Theta_f(5);
        
    f1 = 1 / (sigma_1 * abs(T_deseason(t-1)) * sqrt(2*pi)) ...
        * exp(-1/2 * ((T_deseason(t) - (1 + kappa) * T_deseason(t-1)) ...
        / (sigma_1 * T_deseason(t-1)))^2);

    f2 = 1 / (sigma_2 * sqrt(2*pi)) ...
        * exp(-1/2 * ((T_deseason(t) - (T_deseason(t-1) + mu_2)) ...
        / (sigma_2))^2);
end
%}


