function [ functionLogisticintegral ] = functionLogisticintegral( xlow, xhigh, x0, beta)
% 170630 fcp Integral of logistic function from xlow to xhigh 
%   unnormalized, so not called "cdfLogistic"
%   Notes on fitting 6 MeV spectrum 170629-30

  if (beta==0)
      % beta = 0 is a step function
      % we don't treat this case for now  
  %   nl = length(xlow);
  %   nh = length(xhigh);
  %   for i=1:nl
  %    if (xlow > x0)
  %        if (xhigh > x0)
             functionLogisticintegral = 0.;
  %        else
  %           functionLogisticintegral = xhigh - x0;
  %        end
  %    else
  %        if (xhigh > x0)
  %           functionLogisticintegral = x0 - xlow;
  %        else
  %           functionLogisticintegral = xhigh - xlow;
  %        end
  %    end
  else
      functionLogisticintegral = beta.*log((1 + exp(-(xlow-x0)/beta))./(1 + exp(-(xhigh-x0)/beta)));
  end
end

