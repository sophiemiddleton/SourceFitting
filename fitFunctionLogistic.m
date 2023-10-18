function [ fitFunctionLogistic ] = fitFunctionLogistic(x,par,type)
% fitFunctionLogisticIntegral
%   Return fit function with logistic background and CB lineshape
%   fcp 170630

global ANORM;

% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction 
% 5          CB sigma in channels
% 6          Separation of peaks in channels 
% 7          n (CB power)
% 8          Alpha (CB power law join location)
% 9          Logistic background location parameter
% 10         Logistic background scale parameter
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

if nargin > 2
    % type tells us what component to return; 0 is all coponents
    % 4 is background
    if type > 3
        % background
        fitFunctionLogistic = ANORM(4)*(1-rho)*functionLogistic(x,x0,beta);
        return;
    else
        if type > 0
        % signal type
        fitFunctionLogistic = ANORM(type)*p(type)*rho*functionCB(x, alpha, n, mu(type), sigma);
        return;
        end
    end
end
    % default is add all components
    fitFunctionLogistic = ANORM(4)*(1-rho)*functionLogistic(x,x0,beta) ...
        + ANORM(1)*p(1)*rho*functionCB(x, alpha, n, mu(1), sigma) ...
        + ANORM(2)*p(2)*rho*functionCB(x, alpha, n, mu(2), sigma) ...
        + ANORM(3)*p(3)*rho*functionCB(x, alpha, n, mu(3), sigma);
end

