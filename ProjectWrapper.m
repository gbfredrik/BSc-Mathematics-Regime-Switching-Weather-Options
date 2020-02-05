%% Resetter
% Empties workspace of variables on start
clearvars

%% File loading
% This section does... wanka
sourceDir = "DataSets/";
sourceFiles = ["Stockholm Bromma smhi-opendata 20200130", ...
                "Kiruna Flygplats smhi-opendata 20200130", ...
                "Falsterbo smhi-opendata 20200130"];
%%
if (~exist(sourceDir, 'dir'))
    fprintf("DataSets dir not found!\n\n")
    return
end

for k = 1 : length(sourceFiles)
    %ok loop
    
end

%% Not ok code
if (~exist('data', 'var')) % Only load data once
  [data,txt] = xlsread('labML', 'assetHistory');
  if (size(data,2) == 1)
    S = data(end:-1:1,1);
    dates = datenum(txt(end:-1:3,1));
  else
    S = data(end:-1:1,2);
    dates = datenum(data(end:-1:1,1));
  end
end

%% Data parsing


%%
