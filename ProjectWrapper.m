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
    %A, b, Aeq, beq, lb, ub, nonlcon, options
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

    [Sets(1, k).ML_Theta.VG, Sets(1, k).ML_FVal.VG] = fmincon(@(x) ...
         -log_likelihood_f_VG(Sets(1, k).Deseasoned.Degrees, x(1), x(2), x(3), x(4)), ...
         guess.VG, ...
         [0, -1, -1, 0; 0, -1, 1, 0], [0;0], [], [], [0 0 -inf -inf], [], [], settings.fminconOptions);
    Sets(1, k).ML_FVal.VG = -Sets(1, k).ML_FVal.VG;
     
    [Sets(1, k).ML_Theta.N, Sets(1, k).ML_FVal.N] = fmincon(@(x) ...
        -sum(log(normpdf(Sets(1, k).Deseasoned.Degrees, x(1), x(2)))), ...
        [0 1], [], [], [], [], [], [], [], settings.fminconOptions);
    Sets(1, k).ML_FVal.N = -Sets(1, k).ML_FVal.N;
end
%%
clearvars KS AD GH HYP NIG VG pdf_GH pdf_HYP pdf_NIG pdf_VG cdf_GH cdf_HYP cdf_NIG cdf_VG F

for k = 1 : length(Sets)
    %pdf_series = zeros(length(Sets(1,1).Deseasoned.Degrees),1);
    %sorted_deg = sort(Sets(1,1).Deseasoned.Degrees);
   
    GH = Sets(1,k).ML_Theta.GH;
    HYP = Sets(1,k).ML_Theta.HYP;
    NIG = Sets(1,k).ML_Theta.NIG;
    VG = Sets(1,k).ML_Theta.VG;
    r = -50:0.05:50;
    
    cdf_GH = zeros(length(r), 1);
    cdf_HYP = zeros(length(r), 1);
    cdf_NIG = zeros(length(r), 1);
    cdf_VG = zeros(length(r), 1);
    F = zeros(length(r), 1);
    for i = 1:length(r)
        cdf_GH(i, 1) = integral(@(x) f_GH(x, GH(1), GH(2), GH(3), GH(4), GH(5)), min(r), r(i), 'ArrayValued', true);
        cdf_HYP(i, 1) = integral(@(x) f_GH(x, HYP(1), HYP(2), HYP(3), HYP(4), HYP(5)), min(r), r(i), 'ArrayValued', true);
        cdf_NIG(i, 1) = integral(@(x) f_GH(x, NIG(1), NIG(2), NIG(3), NIG(4), NIG(5)), min(r), r(i), 'ArrayValued', true);
        cdf_VG(i, 1) = integral(@(x) f_VG(x, VG(1), VG(2), VG(3), VG(4)), min(r), r(i), 'ArrayValued', true);
    end

    sorted = sort(Sets(1,k).Deseasoned.Degrees);
    for i = 1 : length(r)
        F(i, 1) = 1 / length(sorted) * sum(sorted <= r(i));
    end

    [KS.h(1), KS.p(1), KS.ks2stat(1)] = kstest2(F, cdf_GH);
    [KS.h(2), KS.p(2), KS.ks2stat(2)] = kstest2(F, cdf_HYP);
    [KS.h(3), KS.p(3), KS.ks2stat(3)] = kstest2(F, cdf_NIG);
    [KS.h(4), KS.p(4), KS.ks2stat(4)] = kstest2(F, cdf_VG);
    [KS.h(5), KS.p(5), KS.ks2stat(5)] = kstest2(F, normcdf(r, Sets(1, k).ML_Theta.N(1), Sets(1, k).ML_Theta.N(2)));
    Sets(1,k).KS = KS;
    
    clearvars KS AD GH HYP NIG VG pdf_GH pdf_HYP pdf_NIG pdf_VG F
end

for k = 1 : length(Sets)
    fprintf('The KS results of set %i are:\n', k);
    disp(Sets(1, k).KS);
end
%%
Sets(1,k).EM = cell(2, 1);
%%
%
clearvars EM Pr1 Pr2 Pr1T Pr2T Theta_f p Q iter_f
close all
model = 1;
useKim = true;

if (model == 1)
    % Initial parameter guesses
    Theta = [0.05, 1, 0.98, 1.5, 1.4, 0.02;
             0.025, 1, 0.98, 10, 1, 0.02;
             0.05, 1, 0.98, 1.5, 1.4, 0.02];

    Theta_Kim = [0.2, 0.1, 0.5, -15, 10, 0.5;
                 0.22, 0.15, 0.8, 15, 1, 0.2;
                 0.1, 1, 0.98, -5, 0.4, 0.02];
elseif (model == 2)
    Theta = [2, 0.5, 1.5, 0.98, -5, 1.5, 0.02;
             0.5, 2, 0.5, 0.98, -3, 1.4, 0.02;
             1, -5, 0.5, 0.98, -20, 1.5, 0.02];

    Theta_Kim = [0.8, -3, 0.5, 0.9, 0, 0.25, 0.1;
                 0.5, -2, 2.5, 0.9, 1, 0.5, 0.1;
                 -0.5, 5, 0.7, 0.9, -1.5, 0.15, 0.1];
end
if (useKim)
   Theta = Theta_Kim; 
end

% Best for model 1:
% p = [0.995 0.005;
%      0.3 0.7];
p = [0.95 0.05;
     0.6 0.4];
% Best for model 2:
% p = [0.995 0.005;
%      0.3 0.7];
% p = [0.995 0.005;
%      0.1 0.9];

for k = 3%1:length(Sets)
    deseasoned = Sets(1,k).Deseasoned.Degrees;% + Sets(1,k).Season_Theta(1) + (1 : length(Sets(1,k).Deseasoned.Degrees))' * Sets(1,k).Season_Theta(2);
    [EM.Pr1, EM.Pr2, EM.Pr1T, EM.Pr2T, EM.Theta_f, EM.p, EM.Q, EM.iter_f] = ...
        EM( deseasoned, ...
            Theta(k,:), ...
            p, ...
            model, ...
            500, ...
            true, ...
            useKim);
    Sets(1,k).EM{model} = EM; clearvars EM
    
    GenerateStatePlot( ...
        Sets(1,k).Deseasoned.Degrees, ...
        Sets(1,k).EM{model}.Pr1T', ...
        Sets(1,k).EM{model}.Pr2T', ...
        0.8);
    
    %figure();
    %plot(Sets(1,k).EM.Q);
 
    fprintf('Finished at iteration: %i, with Q value: %.1f.\n', ...
        Sets(1,k).EM{model}.iter_f, ...
        Sets(1,k).EM{model}.Q(end));
    fprintf('Optimal parameter set is:');
    disp(Sets(1,k).EM{model}.Theta_f(end,:));
    fprintf("P matrix is:\n");
    disp(Sets(1,k).EM{model}.p);
    
    clearvars deseasoned
end

% Clear temp vars
%% Simulate the MRS/HMM
close all
[series, states] = SimulateHMM([0.7895 0.2105], [0.92 0.08; 0.30 0.70], [0.01 0.10], [0.025 0.25], 1, 250, 1)
GenerateStatePlot( ...
    series, ...
    1 - states, ...
    states, ...
    0.5);

%% Simulate and evaluate?
%% Gyamerah
clear G
k = 1;
G.EM = Sets(1,k).EM{1};
par = G.EM.Theta_f(end,:);
par_Levy = Sets(1,k).ML_Theta.GH;

T = 250; %1825
[G.DAT, G.states] = SimulateHMM([G.EM.Pr1T(end), G.EM.Pr2T(end)], ...
    G.EM.p, par(1:3), par(4:6), 2, T, 10, par_Levy', 1);

close all
% GenerateStatePlot( ...
%     G.DAT, ...
%     1 - G.states, ...
%     G.states, ...
%     0.5);

figure(123);
hold on
plot(1:T, G.DAT + seasonFunction(Sets(1,k).Season_Theta, 1:T)')
plot(1:T, seasonFunction(Sets(1,k).Season_Theta, 1:T)')
hold off   

%% Evarest
clear EV
k = 3;
EV.EM = Sets(1,k).EM{2};
par = EV.EM.Theta_f(end,:);

T = 1825; %1825
[EV.DAT, EV.states] = SimulateHMM([EV.EM.Pr1T(end), EV.EM.Pr2T(end)], ...
    EV.EM.p, par(1:4), par(5:7), 3, T, 1000);

close all
% GenerateStatePlot( ...
%     EV.DAT, ...
%     1 - EV.states, ...
%     EV.states, ...
%     0.5);

figure(123);
plot(1:T, EV.DAT)
% 
figure(456)
hold on
plot(1:T, EV.DAT + seasonFunction(Sets(1,k).Season_Theta, 1:T)')
plot(1:T, seasonFunction(Sets(1,k).Season_Theta, 1:T)', 'black')
hold off
xlabel('Days')
ylabel('Daily average temperature')
%%
OOS = Sets(1,3).Clean.Degrees(Sets(1,3).Clean.Time >= Sets(1,3).OutOfSample(1));
OOS = OOS(1:end-1);
figure(10)
hold on
plot(1:T, EV.DAT + seasonFunction(Sets(1,k).Season_Theta, 1:T)', 'blue')
plot(1:T, OOS, 'red')
plot(1:T, seasonFunction(Sets(1,k).Season_Theta, 1:T)', 'black')
hold off
xlabel('Days')
ylabel('Daily average temperature')

%% Calculate HDD/CDD?

OOS_PERIOD = Sets(1,3).Clean.Time(Sets(1,3).Clean.Time >= Sets(1,3).OutOfSample(1));
OOS_PERIOD = OOS_PERIOD(1:end-1);
simulated = EV.DAT + seasonFunction(Sets(1,k).Season_Theta, 1:T)';

JJA = [OOS_PERIOD(month(OOS_PERIOD) == 6); OOS_PERIOD(month(OOS_PERIOD) == 7); OOS_PERIOD(month(OOS_PERIOD) == 8)];
juneDAT = simulated(month(OOS_PERIOD) == 6,:);
julyDAT = simulated(month(OOS_PERIOD) == 7,:);
augDAT = simulated(month(OOS_PERIOD) == 8,:);

T_ref = 18;
for i = 1:1000
    juneCDD{i} = juneDAT(juneDAT(:,i) >= T_ref, i) - T_ref;
    julyCDD{i} = julyDAT(julyDAT(:,i) >= T_ref, i) - T_ref;
    augCDD{i} = augDAT(augDAT(:,i) >= T_ref, i) - T_ref;
end

CDD = 0;
for i = 1:1000
    CDD = CDD + mean(juneCDD{i}) + mean(augCDD{i}) + mean(augCDD{i});
end
CDD_Actual = sum(OOS(OOS >= T_ref) - T_ref)

