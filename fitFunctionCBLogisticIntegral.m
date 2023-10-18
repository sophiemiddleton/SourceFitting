function [ fitFunctionLogisticIntegral ] = fitFunctionLogisticIntegral(par, xl, xh)
% fitFunctionLogisticIntegral
%   Return integral of fit function with logistic background and CB lineshape
%   fcp 170630
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

global ANORM;
global E0;
% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;
global ME;
global MEconstraint;

rho = par(Ifsignal);
alpha = par(IalphaCB);
n = par(InCB);
sigma = par(Isigma);
m = ME*par(Ich613)/E0;   % electron mass in channels
if(MEconstraint)
    mu = [par(Ich613), par(Ich613)-m, par(Ich613)-2.*m];
else
    delta = par(Ich511);
    mu = [par(Ich613), par(Ich613)-delta, par(Ich613)-2.*delta];
end

p = [par(If613), par(Ifescape1), 1-par(If613)-par(Ifescape1)];
x0 = par(IbgLocation);
beta = par(IbgScale);

signal = 0.;
for i = 1:3
    signal = signal + ANORM(i)*p(i)*(functionCBintegral(xl, alpha, n, mu(i), sigma) ...
    - functionCBintegral(xh, alpha, n, mu(i), sigma));
end

fitFunctionLogisticIntegral = rho*signal + ANORM(4)*(1-rho)*functionLogisticintegral(xl,xh,x0,beta);

end

