function [ fitFunctionGComptonNorms ] = fitFunctionGComptonNorms(par, xl, xh)
% fitFunctionCBComptonIntegral
%   Return normalization of fit function components with Compton background and CB lineshape
%   fcp 170926
% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction 
% 5          CB sigma in channels
% 6          Separation of peaks in channels (if not constrained)
% xl, xh     Normalization range

global ME;          % electron mass (any units)
global E0;          % primary photon energy (same units as ME)
global MEconstraint;
% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;

m = ME*par(1)/E0;   % electron mass in channels
sigma = par(Isigma);
if(MEconstraint)
    mu = [par(1), par(1)-m, par(1)-2.*m];
else
    delta = par(Ich511);
    mu = [par(1), par(1)-delta, par(1)-2.*delta];
end
A = zeros(1,4);
for i = 1:3
    A(i) = normcdf(xh, mu(i), sigma) - normcdf(xl, mu(i), sigma);
end

A(4) = functionComptonNormalIntegral(xl,xh,par(1),m,sigma);

fitFunctionGComptonNorms = 1./A;
end

