function [ fitFunctionCBCompton ] = fitFunctionCBCompton(x,par,type)
% fitFunctionLogisticIntegral
%   Return fit function with Compton background and CB lineshape
%   fcp 170915

global ANORM;

% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction  
% 5          CB sigma in channels
% 6          Separation of peaks in channels
% 7          Alpha (CB power law join location)
% 8          n (CB power)
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

p = [par(2), par(3), 1-par(2)-par(3)];

if nargin > 2
    % type tells us what component to return; 0 is all coponents
    % 4 is Compton convoluted with normal with resolution sigma
    if type > 3
        % Compton
        fitFunctionCBCompton = ANORM(4)*(1-rho)*functionComptonNormal(x,par(1),m,sigma);
        return;
    else
        if type > 0
        % signal type
        fitFunctionCBCompton = ANORM(type)*p(type)*rho*functionCB(x, alpha, n, mu(type), sigma);
        return;
        end
    end
end
    % default is add all components
    fitFunctionCBCompton = ANORM(4)*(1-rho)*functionComptonNormal(x,par(1),m,sigma) ...
        + ANORM(1)*p(1)*rho*functionCB(x, alpha, n, mu(1), sigma) ...
        + ANORM(2)*p(2)*rho*functionCB(x, alpha, n, mu(2), sigma) ...
        + ANORM(3)*p(3)*rho*functionCB(x, alpha, n, mu(3), sigma);
end

