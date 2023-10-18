function [ functionomptonintegral ] = functionComptonintegral( xlow, xhigh, x0, beta)
% 170915 fcp Integral of Compton function from xlow to xhigh 
%   unnormalized, so not called "cdfCompton"
%   Notes on fitting 6 MeV spectrum 170629-30, and Compton 170915

  if (beta==0)
      % beta = 0 is a step function
      % we don't treat this case for now  
 
  else
      functionComptonintegral = beta.*log((1 + exp(-(xlow-x0)/beta))./(1 + exp(-(xhigh-x0)/beta)));
  end
end

