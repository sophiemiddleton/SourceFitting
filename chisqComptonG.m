function [ chisq ] = chisqComptonG( par )
% chisqLogistic
%   Return chisq for energy spectrum 
%   Assumes Compton background and Gaussian lineshape
%   fcp 170922

global ANORM;   % fit function normalizations, for each component
global BINCONTENTS;     % Observed bin contents
global EBINS;           % Lower bin energies in channels
global NSAMPLE;
global XLO;     % low end of fit range (channels)
global XHI;     % high end of fit range (channels)

% Normalize each component

ANORM = NSAMPLE*fitFunctionGComptonNorms(par, XLO, XHI);  % Normalization constants

chisq = sum(((BINCONTENTS - fitFunctionGComptonIntegral(par, EBINS, EBINS+1.)).^2)./BINCONTENTS);

end

