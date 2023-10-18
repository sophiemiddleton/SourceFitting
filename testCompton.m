% Script testCompton
% Test Compton functions
% fcp 170915

global ME;          % electron mass (any units)
global E0;          % primary photon energy (same units as ME)
global ANORM;

TESTfunctionCompton = false;
TESTfunctionComptonNormal = false;
TESTfunctionComptonNormalIntegral = false;
TESTfitFunctionCBCompton = true;
TESTfitFunctionCBComptonIntegral = false;

if(TESTfunctionCompton)
    m = 0.511;
    x0 = 6.13;
%    x0 = 0.662;
    % x is KE of Compton electron
    % scattered photon energy is x0-x
    xmax = x0*(1 - 1/(1 + 2*x0/m));
    xmax = x0;
    xinterval = [0,xmax];
   
    channelTEST = false;
    if(channelTEST)    
    %begin test
        % paramter array
    % parameters:
% 1          Energy of primary photon, in channels
par(1) = 194.;
% 2          Fraction of signal in primary photon
par(2) = 0.64;
% 3          Fraction of signal in first escape peak
par(3) = 0.21;
% 4          Signal fraction  
par(4) = 0.52;
par(4) = 0.;
% 5          Alpha (CB power law join location)
par(5) = 0.55;
% 6          n (CB power)
par(6) = 3.8;
% 7          CB sigma in channels
par(7) = 10.;
% 8          Separation of peaks in channels
par(8) = 7.9;
% xl, xh     Integration range


    ME = 0.511;
    E0 = 6.13;
    m = ME*par(1)/E0;   % electron mass in channels
    sigma = par(7);
    x0 = par(1);
    xlo = 149.;
    xhi = 230.;
    NSAMPLE = 1513526;
    % Normalize each component
    ANORM = NSAMPLE*fitFunctionCBComptonNorms(par, xlo, xhi);  % Normalization constants
    xinterval = [xlo,xhi];
    type = 4;
    end
    %end test
    
    
    
    f = @(x)functionCompton(x,x0,m);
    figure;
    fplot(f, xinterval)
    title('TESTfunctionCompton');
%    ylim([0 27]);
%    x = -1:.05:5;
%    xt = x';
%    y=functionLogistic(xt,beta);
end

if(TESTfunctionComptonNormal)
    m = 0.511;
    x0 = 6.13;
    sigma = .5;
%    x0 = 0.662;
    % x is KE of Compton electron
    % scattered photon energy is x0-x
    xmax = x0*(1 - 1/(1 + 2*x0/m));
    xinterval = [0,x0];
    
    
    channelTEST = true;
    if(channelTEST)    
    %begin test
        % paramter array
    % parameters:
% 1          Energy of primary photon, in channels
par(1) = 194.;
%par(1) = 6.13;
% 2          Fraction of signal in primary photon
par(2) = 0.64;
% 3          Fraction of signal in first escape peak
par(3) = 0.21;
% 4          Signal fraction  
par(4) = 0.52;
par(4) = 0.;
% 5          Alpha (CB power law join location)
par(5) = 0.55;
% 6          n (CB power)
par(6) = 3.8;
% 7          CB sigma in channels
par(7) = 10.;
% 8          Separation of peaks in channels
par(8) = 7.9;
% xl, xh     Integration range


    ME = 0.511;
    E0 = 6.13;
    m = ME*par(1)/E0;   % electron mass in channels
    sigma = par(7);
    x0 = par(1);
    xlo = 149.;
    xhi = 230.;
    xlo = 130;
    xhi = 230; 
    
    NSAMPLE = 1513526;
    % Normalize each component
    %ANORM = NSAMPLE*fitFunctionCBComptonNorms(par, xlo, xhi);  % Normalization constants
    xinterval = [xlo,xhi];
    type = 4;
    end
    %end test
    
    f = @(x)functionComptonNormal( x, x0, m, sigma);
    figure;
    fplot(f, xinterval)
  %  ylim([0 8]);
    title('TESTfunctionComptonNormal');
end
   
if(TESTfunctionComptonNormalIntegral)
    m = 0.511;
    x0 = 6.13;
    sigma = 0.5;
    xlo = 2.;
    xinterval = [xlo,1.3*x0];
    f = @(x)functionComptonNormalIntegral(xlo,x,x0,m,sigma);
    figure;
    fplot(f, xinterval)
    title('TESTfunctionComptonNormalIntegral');
end

if(TESTfitFunctionCBCompton)
    % paramter array
    % parameters:
% 1          Energy of primary photon, in channels
par(1) = 194.;
% 2          Fraction of signal in primary photon
par(2) = 0.64;
% 3          Fraction of signal in first escape peak
par(3) = 0.21;
% 4          Signal fraction  
par(4) = 0.52;
%par(4) = 0.;
% 5          Alpha (CB power law join location)
par(5) = 0.55;
% 6          n (CB power)
par(6) = 3.8;
% 7          CB sigma in channels
par(7) = 10.;
% 8          Separation of peaks in channels
par(8) = 7.9;
% xl, xh     Integration range


    ME = 0.511;
    E0 = 6.13;
    m = ME*par(1)/E0;   % electron mass in channels
    sigma = par(7);
    xlo = 149.;
    xhi = 230.;
    NSAMPLE = 1513526;
    % Normalize each component
    ANORM = NSAMPLE*fitFunctionCBComptonNorms(par, xlo, xhi);  % Normalization constants
    xinterval = [xlo,xhi];
    %
    type = 0;
    f = @(x)fitFunctionCBCompton(x,par,type);
    %f = @(x)fitFunctionCBCompton(x,par);
    figure;
    fplot(f, xinterval)
    title('TESTfitFunctionCBCompton');
end

if(TESTfitFunctionCBComptonIntegral)
    % paramter array
    % parameters:
% 1          Energy of primary photon, in channels
par(1) = 194.;
% 2          Fraction of signal in primary photon
par(2) = 0.64;
% 3          Fraction of signal in first escape peak
par(3) = 0.21;
% 4          Signal fraction  
par(4) = 0.52;
par(4) = 0.;
% 5          Alpha (CB power law join location)
par(5) = 0.55;
% 6          n (CB power)
par(6) = 3.8;
% 7          CB sigma in channels
par(7) = 10.;
% 8          Separation of peaks in channels
par(8) = 7.9;
% xl, xh     Integration range
    ME = 0.511;
    E0 = 6.13;
    m = ME*par(1)/E0;   % electron mass in channels
    sigma = par(7);
    xlo = 149.;
    xhi = 230.;
    NSAMPLE = 1513526;
    % Normalize each component
    ANORM = NSAMPLE*fitFunctionCBComptonNorms(par, xlo, xhi);  % Normalization constants

    xinterval = [xlo,xhi];
    f = @(x)fitFunctionCBComptonIntegral(par, xlo, x);
    figure;
    fplot(f, xinterval)
    title('TESTfitFunctionCBComptonIntegral');
end