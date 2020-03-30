classdef Settings
    % Class definition of Settings class
    %   Matlab doesn't support constant variables in any pretty format,
    %   so here's a workaround by 
   properties (Constant)
      N = 10000
      %D = 1/NamedConst.R
      avgType = 'MinMax' % Set average type DAT series, toggle 'MinMax' or ' Mean'
      
   end
end