function [ functionLogistic ] = functionLogistic( x, x0, beta )
% 170630 fcp logistic function
%   unnormailzed, so not called "pdfLogistic"
%   Notes on fitting 6 MeV spectrum 170629-30
 if(beta == 0)
     % beta = 0 is a step function
     ndo = length(x);
     for i=1:ndo
     if (x(i) > x0) 
         functionLogistic(i) = 0.;
     else
         functionLogistic(i) = 1.;
     end
     end
 else
     functionLogistic = 1./(1+exp((x-x0)/beta));
 end
end

