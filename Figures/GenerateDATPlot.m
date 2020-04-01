function [status] = GenerateDATPlot(Set, seasonFunction, X, ...
    showFigures, saveFigures, showTref, showLinTrend)

try
    if showFigures
        f = figure();
    else
        f = figure('visible', 'off');
    end
    
    hold on
    title(Set.ShortName)
    plot(Set.CleanSet.Time(end-3650:end,:), transpose(Set.CleanSet.Degrees(end-3650:end,:)), 'b')
    plot(seasonFunction(X(1), X(2), X(3), X(4), 0:3650), 'r', 'LineWidth', 2)
    legendLabels = ["DAT", "Seasonal"];

    if showTref
        plot(18 .* ones(3650, 1), 'g', 'LineWidth', 2)
        legendLabels = [legendLabels, "T_r_e_f"];
    end

    if showLinTrend
        plot(seasonFunction(X(1), X(2), 0, 0, 0:3650), 'c', 'LineWidth', 2)
        legendLabels = [legendLabels, "Linear trend"];
    end
    legend(legendLabels, 'Location', 'EastOutside')
    hold off
    
    if saveFigures
        %print(f, sprintf('Figures/DAT/%s DAT plot %s', data.ShortName, datestr(now)), '-dpng');
        print(f, sprintf('Figures/DAT/%s DAT Plot', Set.ShortName), '-dpng');
    end
    status = 1;
catch
    fprintf('An error occured in the DAT plot function.\n')
    status = 0;
end

end
