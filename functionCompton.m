function [ functionCompton ] = functionCompton( x, x0, m )
% 170915 fcp Compton function
%   unnormalized, so not called "pdfCompton"
%   Notes on fitting 6 MeV spectrum 170629-30, plus Compton 170915
%   Return value is proportional to the density function for 
%   the KE of the Compton electron
%   x is the KE of the Compton electron
%   x0 is the incident photon energy
%   m is normally electron mass
%   Compute maximum possible x:
    xmax = x0*(1 - 1/(1 + 2*x0/m));
    belowmax = x < xmax;
    r = x0./(x0-x);
    cos = 1 - (m/x0)*x./(x0-x);
%   functionCompton(belowMax) = r + 1./r - 1. + cos.^2;
%    functionCompton(~belowMax) = 0.;
 functionCompton = belowmax.*(r + 1./r - 1. + cos.^2) + (~belowmax)*0;
end

