function y = ghpdf(x,Param,m,sigma)
%GHPDF Generalized Hyperbolic probability density function (pdf).
%   Y = GHPDF(X,PARAM,M,SIGMA) Returns the GH pdf
%   with parameters: PARAM = [LAMBDA, ZETA, RHO], mean M, and standard
%   deviation SIGMA at the values in X. 2nd parametrization in Prauss
%   (1999).
%
%   LAMBDA in IR intensifies the skewness (when RHO <> 0)
%   ZETA > 0 is a shape parameter
%   -1 < RHO < 1 is the skewness coefficient
%
%   When RHO > 0, the GHD becmes skewed to the right, and when RHO < 0, it
%   becomes skewed to the left.
%
%   The size of Y is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Default values for PARAM and M are [0, Inf, 0], 0, and 1 respectively.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Copyright (c) 22 November 2014 by Ahmed BenSaïda           %
%                 LaREMFiQ Laboratory, IHEC Sousse - Tunisia             %
%                       Email: ahmedbensaida@yahoo.com                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1, 
    error('Requires at least one input argument.');
end
if nargin < 2 || isempty(Param)
    Param = [0, 500, 0]; % These are the default values for te normal distribution.
end
if size(Param,2) ~= 3
    error('Input argument ''Param'' must be a 3-columns parameter.')
end
if nargin < 3 || isempty(m);
    m = 0;
end
if nargin < 4 || isempty(sigma),
    sigma = 1;
end
% Extract distribution parameters from PARAM.
Lambda = Param(:,1);
Zeta   = Param(:,2);
Rho    = Param(:,3);
[errorcode, x, Lambda, Zeta, Rho, m, sigma] = paramchck(6,x,Lambda,Zeta,Rho,m,sigma);
if errorcode > 0
    error('Requires non-scalar arguments to match in size.');
end
%   Initialize Y to zero.
if isa(x,'single') || isa(Zeta,'single') || isa(Rho,'single') || ...
        isa(Lambda,'single') || isa(m,'single') || isa(sigma,'single')
    y = zeros(size(x),'single');
else
    y = zeros(size(x));
end
% Return NaN for invalid parameters.
y(Zeta <= 0 | abs(Rho) >= 1 | sigma <= 0 | isnan(x) | isnan(m) | ...
    isnan(sigma) | isnan(Lambda) | isnan(Zeta) | isnan(Rho)) = NaN;
h   =   find(Zeta > 0 & abs(Rho) < 1 & sigma > 0);
if any(h)
    Delta = 1 ./ sqrt(besselk(Lambda(h) + 1, Zeta(h)) ./ (Zeta(h) .* ...
        besselk(Lambda(h), Zeta(h))) + (Rho(h).^2 ./ (1 - Rho(h).^2)) .* ...
        (besselk(Lambda(h) + 2, Zeta(h)) ./ besselk(Lambda(h), Zeta(h)) - ...
        (besselk(Lambda(h) + 1, Zeta(h)) ./ besselk(Lambda(h), Zeta(h))).^2));
    
    mu   = m(h) -  sigma(h) .* Delta .* Rho(h) ./ sqrt(1 - Rho(h).^2) .* ...
        besselk(Lambda(h) + 1, Zeta(h)) ./ besselk(Lambda(h), Zeta(h));
    
    xx   = ((x(h) - mu) ./ (sigma(h) .* Delta));
    
    y(h) =  (1 - Rho(h).^2).^(Lambda(h)/2 - 0.25) .* sqrt(Zeta(h)) ./ ...
        (sqrt(2*pi) * sigma(h) .* Delta .* besselk(Lambda(h), Zeta(h))) .* ...
        (1 + xx.^2).^ (Lambda(h)/2 - 0.25) .* besselk(Lambda(h) - 1/2, ...
        Zeta(h) ./ sqrt(1 - Rho(h).^2) .* sqrt(1 + xx.^2)) .* ...
        exp(Zeta(h) .* Rho(h) ./ sqrt(1 - Rho(h).^2) .* xx);
end