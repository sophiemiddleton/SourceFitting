function [ fitFunctionCBLogisticNorms ] = fitFunctionCBLogisticNorms(par, xl, xh)
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
% 9          Logistic background location parameter
% 10         Logistic background scale parameter
% xl, xh     Normalization range

global ME;          % electron mass (any units)
global E0;          % primary photon energy (same units as ME)
global MEconstraint;
% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;

m = ME*par(Ich613)/E0;   % electron mass in channels
rho = par(Ifsignal);
sigma = par(Isigma);
p = [par(If613), par(Ifescape1), 1-par(If613)-par(Ifescape1)];
if(MEconstraint)
    mu = [par(Ich613), par(Ich613)-m, par(Ich613)-2.*m];
else
    delta = par(Ich511);
    mu = [par(Ich613), par(Ich613)-delta, par(Ich613)-2.*delta];
end

alpha = par(IalphaCB);
n = par(InCB);
p = [par(If613), par(Ifescape1), 1-par(If613)-par(Ifescape1)];
x0 = par(IbgLocation);
beta = par(IbgScale);

A = zeros(1,4);
for i = 1:3
    A(i) = functionCBintegral(xl, alpha, n, mu(i), sigma) ...
    - functionCBintegral(xh, alpha, n, mu(i), sigma);
end

A(4) = functionLogisticintegral(xl,xh,x0,beta);

fitFunctionCBLogisticNorms = 1./A;
end

