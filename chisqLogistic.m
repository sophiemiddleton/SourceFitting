function [ chisq ] = chisqLogistic( par )
% chisqLogistic
%   Return chisq for energy spectrum 
%   Assumes logistic background and Crystal Ball lineshape
%   fcp 170630

global ANORM;   % fit function normalization
global BINCONTENTS;     % Observed bin contents
global EBINS;           % Lower bin energies in channels
global NSAMPLE;
global XLO;     % low end of fit range (channels)
global XHI;     % high end of fit range (channels)


ANORM = NSAMPLE*fitFunctionCBLogisticNorms(par, XLO, XHI);  % Normalization constants

chisq = sum(((BINCONTENTS - fitFunctionCBLogisticIntegral(par, EBINS, EBINS+1.)).^2)./BINCONTENTS);

end

