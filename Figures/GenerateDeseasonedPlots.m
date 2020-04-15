function [status] = GenerateDeseasonedPlots(Set, showFigures, saveFigures)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% try
    if showFigures
        f1 = figure();
    else
        f1 = figure('visible', 'off');
    end
    plot(Set.Deseasoned.Time, Set.Deseasoned.Degrees)
    ylabel('Degrees Celsius')
    box off

    if showFigures
        f2 = figure();
    else
        f2 = figure('visible', 'off');
    end
    qqplot(Set.Deseasoned.Degrees)
    title('')
    if saveFigures
        %print(f, sprintf('Figures/DAT/%s DAT plot %s', data.ShortName, datestr(now)), '-dpng');
        if verLessThan('matlab', '9.8.0')
            print(f1, sprintf('Figures/Deseasoned/%s Deseason', Set.ShortName), '-dpng');
            print(f2, sprintf('Figures/Deseasoned/%s QQ', Set.ShortName), '-dpng');
        else
            exportgraphics(f1, sprintf('Figures/Deseasoned/%s Deseason%s', Set.ShortName, '.png'));
            exportgraphics(f2, sprintf('Figures/Deseasoned/%s QQ%s', Set.ShortName, '.png'));
        end
    end
    status = 1;
% catch
%      fprintf('An error occured in the Deseasoned plot function.\n')
%      status = 0;
% end

end

