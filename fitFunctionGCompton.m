function [ fitFunctionGCompton ] = fitFunctionGCompton(x,par,type)
% fitFunctionGCompton
%   Return fit function with Compton background and Gaussian lineshape
%   fcp 170915

global ANORM;
global MEconstraint;

% parameters:
% 1          Energy of primary photon, in channels
% 2          Fraction of signal in primary photon
% 3          Fraction of signal in first escape peak
% 4          Signal fraction  
% 5          CB sigma in channels
% 6          Separation of peaks in channels unless comstrained

% need electron mass in channels
global ME;          % electron mass (any units)
global E0;          % primary photon energy (same units as ME)

m = ME*par(1)/E0;   % electron mass in channels
rho = par(4);
sigma = par(5);
p = [par(2), par(3), 1-par(2)-par(3)];
if(MEconstraint)
    mu = [par(1), par(1)-m, par(1)-2.*m];
else
    delta = par(6);
    mu = [par(1), par(1)-delta, par(1)-2.*delta];
end


if nargin > 2
    % type tells us what component to return; 0 is all coponents
    % 4 is Compton convoluted with normal with resolution sigma
    if type > 3
        % Compton
        fitFunctionGCompton = ANORM(4)*(1-rho)*functionComptonNormal(x,par(1),m,sigma);
        return;
    else
        if type > 0
        % signal type
        fitFunctionGCompton = ANORM(type)*p(type)*rho*normpdf(x, mu(type), sigma);
        return;
        end
    end
end
    % default is add all components
    fitFunctionGCompton = ANORM(4)*(1-rho)*functionComptonNormal(x,par(1),m,sigma) ...
        + ANORM(1)*p(1)*rho*normpdf(x, mu(1), sigma) ...
        + ANORM(2)*p(2)*rho*normpdf(x, mu(2), sigma) ...
        + ANORM(3)*p(3)*rho*normpdf(x, mu(3), sigma);
end

