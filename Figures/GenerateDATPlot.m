function [status] = GenerateDATPlot(Set, seasonFunction, X, ...
    showFigures, saveFigures, showSeason, showTref, showLinTrend, setPeriod)

try
    if showFigures
        f = figure();
    else
        f = figure('visible', 'off');
    end
    
    if setPeriod == "In"
        period = Set.InSample;
    elseif setPeriod == "InOut"
        period = [Set.InSample, Set.OutOfSample];
    end
    
    hold on
    title(Set.ShortName)
    plot(Set.Clean.Time(period), transpose(Set.Clean.Degrees(period)), 'b')
    legendLabels = "DAT";
    ylabel('Degrees Celsius')
    
    if showSeason
        plot(seasonFunction(X, 0:length(period)-1), 'r', 'LineWidth', 2)
        legendLabels = [legendLabels, "Seasonal"];
    end
    
    
    if showTref
        plot(18 .* ones(length(period), 1), 'g', 'LineWidth', 2)
        legendLabels = [legendLabels, "T_r_e_f"];
    end

    if showLinTrend
        plot(seasonFunction([X(1), X(2), 0, 0], 0:length(period)-1), 'c', 'LineWidth', 2)
        legendLabels = [legendLabels, "Linear trend"];
    end
    legend(legendLabels, 'Location', 'EastOutside')
    hold off
    
    if saveFigures
        %print(f, sprintf('Figures/DAT/%s DAT plot %s', data.ShortName, datestr(now)), '-dpng');
        if verLessThan('matlab', '9.8.0')
            print(f, sprintf('Figures/DAT/%s DAT %s', Set.ShortName, setPeriod), '-dpng');
        else
            exportgraphics(f, sprintf('Figures/DAT/%s DAT %s%s', Set.ShortName, setPeriod, '.png'));
        end
    end
    status = 1;
catch
     fprintf('An error occured in the DAT plot function.\n')
     status = 0;
end

end
