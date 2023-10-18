function [ functionCBintegral ] = functionCBintegral(  x, alpha, n, mu, sigma)
% 170629 fcp Integral of Crystal Ball function 
%   unnormalized, so not called "cdfCB"
%   Notes on fitting 6 MeV spectrum 170629
%   See also crystalball_integral in ROOT 

% If UPPER, the integral is the right tail integral.
% This function may not work if alpha < 0, meaning power law is on right (could be checked).     
% parameters:
% alpha : is non equal to zero, define the # of sigma from which it becomes a power-law function (from mean-alpha*sigma)
% if alpha > 0, power law is on left, and int_x^infty is computed.
% if alpha < 0, power law is on right, and int_-infty^x is computed.
% n > 1 : is integrer, is the power of the low  tail
% add a value xmin for cases when n <=1 the integral diverges 
  if (alpha==0)
    fprintf('functionCBintegral - CrystalBall function not defined at alpha=0\n');
    return;
  end   
  if (n<=0)   
    fprintf('functionCBintegral - n<=0 not allowed\n');
    return;
  end  
%  optional UPPER argument to decide which way to integrate
% if nargin < 6
%    UPPER = 0;
% end
 z = sign(alpha)*(x-mu)/sigma;
 ndo = length(z);
 functionCBintegral = zeros(1,ndo);
 if (sigma <= 0.)   
     return 
 end

  useLog = (n == 1.0); 
  abs_alpha = abs(alpha);
   
  
  sqrtpiover2 = sqrt(pi/2.);
  sqrt2pi = 2.*sqrtpiover2; 
  oneoversqrt2 = 1./sqrt(2.);
  
  for ido=1:ndo
      if (z(ido) <= -abs_alpha)
        A = (n/abs_alpha)^n * exp(-0.5 * alpha^2);
        B = n/abs_alpha - abs_alpha;
           if (~useLog) 
              C = (n/abs_alpha) * (1./(n-1)) * exp(-0.5*alpha^2);
              intpow  = C - A /(n-1.) * (B-z(ido))^(-n+1);
           else
               % for n=1 the primitive of 1/x is log(x)
               intpow = A * log( (B - z(ido)) * abs_alpha/n );
           end
            intgaus = sqrtpiover2*(1. + erf(abs_alpha*oneoversqrt2));
      else
            intgaus = sqrt2pi*normcdf(z(ido), 0, 1, 'upper');
            intpow  =  0;  
      end
      functionCBintegral(ido) = sigma * (intgaus + intpow);
  end
end

