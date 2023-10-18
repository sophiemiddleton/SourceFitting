function [ functionComptonNormal ] = functionComptonNormal( x, x0, m, sigma)
% 170920 fcp Compton function convoluted with a normal
%   unnormalized, so not called "pdfComptonNormal"
%   Notes on fitting 6 MeV spectrum 170629-30, plus Compton 170915-21
%   Return value is proportional to the density function for 
%   the KE of the Compton electron convoluted with a normal
%   x is the KE of the Compton electron as observed
%   x0 is the incident photon energy
%   m is normally electron mass
%   sigma is standard deviation of the normal
%   Compute maximum possible electron KE:
    xmax = x0*(1 - 1/(1 + 2*x0/m));
%    r = x0./(x0-x);
%    cos = 1 - (m/x0)*x./(x0-x);
 n = length(x);
 A = zeros(1,n);
 for in = 1:n
 fun = @(y)normpdf(y-x(in), 0, sigma).*functionCompton(y, x0, m);
 A(in) =  integral(fun, 0, xmax);
 end
 functionComptonNormal = A;
end

