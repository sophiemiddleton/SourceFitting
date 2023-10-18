function [ functionComptonNormalIntegral ] = functionComptonNormalIntegral( xlo, xhi, x0, m, sigma)
% 170921 fcp Integral of Compton function convoluted with a normal
%   unnormalized, so not called "cdfComptonNormal"
%   Notes on fitting 6 MeV spectrum 170629-30, plus Compton 170915-21
%   Return value is proportional to the density function for 
%   the KE of the Compton electron convoluted with a normal
%   x is the KE of the Compton electron as observed
%   x0 is the incident photon energy
%   m is normally electron mass
%   sigma is standard deviation of the normal
%   Compute maximum possible electron KE:
%    xmax = x0*(1 - 1/(1 + 2*x0/m));
%    r = x0./(x0-x);
%    cos = 1 - (m/x0)*x./(x0-x);
%   functionCompton(belowMax) = r + 1./r - 1. + cos.^2;
%    functionCompton(~belowMax) = 0.;
 fun = @(x)functionComptonNormal(x, x0, m, sigma);
 N = length(xlo);
 M = length(xhi);
 if (N > M)
     xl = xlo;
     xh(1:N) = xhi;
     T = N;
 elseif (M > N)
     xh = xhi;
     xl(1:M) = xlo;
     T = M;
 else
     xl = xlo;
     xh = xhi;
     T = N;
 end
 A = zeros(1,T);
 for n=1:T 
    A(n) = integral(fun, xl(n), xh(n));
 end
 functionComptonNormalIntegral = A;
end

