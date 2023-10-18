function [ fitFunctionGComptonIntegral ] = fitFunctionGComptonIntegral(par, xl, xh)
% fitFunctionCBComptonIntegral
%   Return integral of fit function with Compton background and Gaussian lineshape
%   fcp 170929
% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction 
% 5          CB sigma in channels
% 6          Separation of peaks in channels unless constrained
% xl, xh     Integration range

global ANORM;
global ME;          % electron mass (any units)
global MEconstraint;
global E0;          % primary photon energy (same units as ME)
% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;

m = ME*par(Ich613)/E0;   % electron mass in channels
rho = par(Ifsignal);
sigma = par(Isigma);
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

signal = 0.;
for i = 1:3
    signal = signal + ANORM(i)*p(i)*(normcdf(xh, mu(i), sigma) ...
    - normcdf(xl, mu(i), sigma));
end

fitFunctionGComptonIntegral = rho*signal +(1-rho)*ANORM(4)*functionComptonNormalIntegral(xl,xh,par(1),m,sigma);

end

