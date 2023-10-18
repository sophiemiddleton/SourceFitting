function [ fitFunctionCBComptonNorms ] = fitFunctionCBComptonNorms(par, xl, xh)
% fitFunctionCBComptonIntegral
%   Return normalization of fit function components with Compton background and CB lineshape
%   fcp 170926
% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction  
% 5          CB sigma in channels
% 6          Separation of peaks in channels
% 7          Alpha (CB power law join location)
% 8          n (CB power)
% xl, xh     Normalization range

% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;
global MEconstraint;
global ME;          % electron mass (any units)
global E0;          % primary photon energy (same units as ME)

rho = par(Ifsignal);
alpha = par(IalphaCB);
n = par(InCB);
sigma = par(Isigma);
m = ME*par(Ich613)/E0;   % electron mass in channels
if(MEconstraint)
    mu = [par(1), par(1)-m, par(1)-2.*m];
else
    delta = par(Ich511);
    mu = [par(1), par(1)-delta, par(1)-2.*delta];
end

A = zeros(1,4);
for i = 1:3
    A(i) = functionCBintegral(xl, alpha, n, mu(i), sigma) ...
    - functionCBintegral(xh, alpha, n, mu(i), sigma);
end

A(4) = functionComptonNormalIntegral(xl,xh,par(1),m,sigma);

fitFunctionCBComptonNorms = 1./A;
end

