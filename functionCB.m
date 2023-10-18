function [ functionCB ] = functionCB( x, alpha, n, mu, sigma )
% 170629 fcp Crystal Ball function
%   unnormalized, so not called "pdfCB"
%   Notes on fitting 6 MeV spectrum 170629
%   See also crystalball_function in ROOT

 if (sigma < 0.)   
     functionCB = 0.;
     return 
 end

 z = sign(alpha)*(x - mu)/sigma; 
 abs_alpha = abs(alpha);
 
 test = z  > - abs_alpha;
 functionCB = test.*exp(- 0.5 * z.^2) + (~test).*(exp(-0.5*alpha^2)*(1-alpha^2/n-z*abs_alpha/n).^(-n));
 
end

