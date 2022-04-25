classdef Settings
    % Class definition of Settings class
    %   Matlab doesn't support constant variables in any pretty format,
    %   so here's a workaround by using a class with constant properties.
    %   Contains: 
    %       N = Number of simulations
    %       avgType = Method used for DAT calculation
    %       fminconOptions = Options for optimization solver
    properties (Constant)
        N = 10000
        %D = 1/NamedConst.R
        avgType = 'MinMax' % Set average type DAT series, toggle 'MinMax' or ' Mean'
        fminconOptions = optimoptions('fmincon', 'Display', 'iter');%, 'StepTolerance', 1e-10, 'FunctionTolerance', 1e-10, 'ConstraintTolerance', 1e-10)
        rndSeed = 'Default'
    end
end
