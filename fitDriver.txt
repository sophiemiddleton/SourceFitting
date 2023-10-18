% Fit source calibration spectrum for mu2e
% Fit it, and determine precision as function of 
% sample size.
% fcp 170616
% fcp mod 170630 chisq fit with logistic background, crystal ball lineshape
% fcp mod 170915 implement swtich for different fit function options
% fcp mod 170922 add option to fit with Compton "background"
% fcp mod 171001 reorder parameters so that sigma is fifth and delta is
% sixth
% fcp mod 171006 Add a pedestal and possibility to constrain ME
%   
%           Declare globals
global ANORM
global BINCONTENTS;
global E0;
E0 = 6.13;
global EBINS;
global EDATA;
global ME;
ME = 0.511;
global NBINS;
global NSAMPLE;
global XHI;
global XLO;
global MEconstraint;
% parameter indices
global Ich613 If613 Ifescape1 Ifsignal Isigma Ich511 IalphaCB InCB IbgLocation IbgScale;

Compton = false;
logistic = ~Compton;
CB = true;
Gauss = ~CB;
MEconstraint = true;

E1 = E0 - ME;
E2 = E1 - ME;

% Read spectrum
delimiterIn = '\t';
headerlinesIn = 0;
spectrum = importdata('dtgen170616CsISIPM.dat',delimiterIn);

x = spectrum(:,1);
y = spectrum(:,2);
nch = length(x);

plot(x,y);
xlim([0,20000]);
figure;

pedFull = 731-1;          % Peak of dtgen170616CsISIPM.dat spectrum occurs here
%pedFull = 1000.;
nbits = 11;
maxCount = 2^nbits-1;
fullScale = 60;
leastCount = fullScale/maxCount;

Epeak = 6.13-0.511;
peakChannel = 14000;
MeVperChannel = Epeak/peakChannel;
nCombine = round(leastCount/MeVperChannel);
nmax = floor(nch/nCombine)*nCombine;
nc = nmax/nCombine;

xCombine = (1:nc)';
yCombine = zeros(nc,1);

for i = 1:nmax
    j = floor((i-1)/nCombine) + 1;
    yCombine(j) = yCombine(j) + y(i);
end

plot(xCombine, yCombine);
xlim([0,300]);
%figure;

pedCombine = pedFull/nCombine;          % Pedestal in 11 bit ADC channels
xCombine = xCombine - pedCombine;       % Pedestal subtraction

% Do a maximum likelihood fit to the data

nsample = 0;

%binstart = 150;        % pre 170928
binstart = 150;
XLO = binstart-1. - pedCombine;      % Bins are counted starting at channel 1, for interval 0 to 1
                                     % XLO and XHI are pedestal-subtracted channel space 
binend   = 230;
XHI = binend - pedCombine;
NBINS = binend - binstart + 1;

maxContents = 0;
for ib = binstart:binend
    % Find bin with maximum counts
    if (yCombine(ib) > maxContents)
        binMax = ib;
        maxContents = yCombine(ib);
    end
    % Set up a vector of data at bin centers to make histogram later
    if (yCombine(ib)>0)
        EDATA(nsample+1:nsample+yCombine(ib)) = xCombine(ib)-.5;
        nsample = nsample + yCombine(ib);
    end
end

for ib = binMax:binend
%    fprintf('%d %d\n',ib,yCombine(ib));
end

EBINS = binstart-1:binend-1;        % Bin energies in channels, shifted for bin 1 to start at 0
EBINS = EBINS - pedCombine;         % Pedestal subtraction
BINCONTENTS = yCombine(binstart:binend)';


%   Set up parameter vector with starting values
par0(1) = 184.3; %E0/(nCombine*MeVperChannel);    % Energy of primary photon, in channels
parname(1) = string('ch613');
Ich613 = 1;
par0(2) = 0.83;         % Fraction of signal in primary photon
parname(2) = 'f613';
If613 = 2;
par0(3) = 0.12;         % Fraction of signal in first escape peak
parname(3) = 'fescape1';
Ifescape1 = 3;
par0(4) = 0.58; %0.5;          % Signal fraction  
parname(4) = 'fsignal';
Ifsignal = 4;
par0(5) = 0.05*E0/(nCombine*MeVperChannel);    % normal/CB sigma in channels
par0(5) = 9.7;
parname(5) = 'sigma';
Isigma = 5;
npar = 5;
if(~MEconstraint)
par0(6) = .511/(nCombine*MeVperChannel);       % Separation of peaks in channels
par0(6) = 15.;
parname(6) = 'ch.511';
Ich511 = 6;
npar = 6;
end
if(CB)
    npar = npar + 1;
par0(npar) = 0.5; %20;            % Alpha (CB power law join location)
parname(npar) = 'alphaCB';
IalphaCB = npar;
    npar = npar + 1;
par0(npar) = 10.;           % n (CB power)
parname(npar) = 'nCB';
InCB = npar;
if(logistic)
    npar = npar + 1;
par0(npar) = 0; %100;            % Logistic background location parameter
parname(npar) = 'bgLocation';
IbgLocation = npar;
    npar = npar + 1;
par0(npar) = 40.; %50         % Logistic background scale parameter
parname(npar) = 'bgScale';
IbgScale = npar;
end
elseif(Gauss)
if(logistic)
    npar = npar + 1;
par0(npar) = 0; %100;            % Logistic background location parameter
parname(npar) = 'bgLocation';
IbgLocation = npar;
    npar = npar + 1;
par0(npar) = 40.; %50         % Logistic background scale parameter
parname(npar) = 'bgScale';
IbgScale = npar;
end
end

ndof = binend-binstart+1-npar;

NSAMPLE = nsample;
    
if(Compton)     % Compton background function
    if(CB)      % Crystal Ball signal function
[parhat, chisq] = fminsearch(@chisqCompton, par0);
[parhatu, chisqu,exitflag,output,grad,hessian] = fminunc(@chisqCompton, parhat);
sigma = parhatu(5);
ch511 = parhatu(6);
sigmaRelative = sigma/ch613;
    elseif(Gauss)   % Normal signal function
[parhat, chisq] = fminsearch(@chisqComptonG, par0);
[parhatu, chisqu,exitflag,output,grad,hessian] = fminunc(@chisqComptonG, parhat);
sigma = parhatu(5);
ch511 = parhatu(6);
    end
else        % Logistic background function
            % Only signal option is Crystal Ball
[parhat, chisq] = fminsearch(@chisqLogistic, par0);

[parhatu, chisqu,exitflag,output,grad,hessian] = fminunc(@chisqLogistic, parhat);
residual = BINCONTENTS - fitFunctionCBLogisticIntegral(parhatu, EBINS, EBINS+1.);
residualNormalized = (BINCONTENTS - fitFunctionCBLogisticIntegral(parhatu, EBINS, EBINS+1.))./sqrt(BINCONTENTS);
chisqA = sum(residualNormalized.^2);

end

ch613 = parhatu(1);
sigma = parhatu(5);
ch511 = parhatu(6);
sigmaRelative = sigma/ch613;
% compute implied pedestal
r613511 = 6.13/.511;
predict613 = r613511*ch511;
predictped = ch613 - predict613;

p = chi2cdf(chisqu,ndof,'upper');
Hinv = (0.5.*hessian)^-1;
sigmap = sqrt(diag(Hinv));
smatrix = diag(1./sigmap);
rho = smatrix*Hinv*smatrix;
ped = ch613 - (6.13/.511)*ch511;

table(parname', parhatu', sigmap)

xmin = XLO;
xmax = XHI;
nx = NBINS;
dx = 1;
edges = xmin:dx:xmax;

h=histc(EDATA,edges);
figure;
%y = dx*fittedFunction(edges,parhat);
%plot(edges,y, 'r','LineWidth',3.);

if(Compton)
    if(CB)
ANORM = NSAMPLE*fitFunctionCBComptonNorms(parhatu, XLO, XHI);  % Normalization constants
fplot(@(x)fitFunctionCBCompton(x,parhatu), [xmin xmax],'r','LineWidth',3.);
    elseif(Gauss)
        ANORM = NSAMPLE*fitFunctionGComptonNorms(parhatu, XLO, XHI);  % Normalization constants
        fplot(@(x)fitFunctionGCompton(x,parhatu), [xmin xmax],'r','LineWidth',3.);
    end
else
ANORM = NSAMPLE*fitFunctionCBLogisticNorms(parhatu, XLO, XHI);  % Normalization constant
fplot(@(x)fitFunctionLogistic(x,parhatu), [xmin xmax],'r','LineWidth',3.);
end
set(gca,'FontSize',14);
hold;
b = bar(edges,h,'histc');

if(Compton)
    if(CB)
fplot(@(x)fitFunctionCBCompton(x,parhatu,1), [xmin xmax],'b--','LineWidth',1.5);
fplot(@(x)fitFunctionCBCompton(x,parhatu,2), [xmin xmax],'g--','LineWidth',1.5);
fplot(@(x)fitFunctionCBCompton(x,parhatu,3), [xmin xmax],'c--','LineWidth',1.5);
fplot(@(x)fitFunctionCBCompton(x,parhatu,4), [xmin xmax],'y--','LineWidth',1.5);
    elseif(Gauss)
fplot(@(x)fitFunctionGCompton(x,parhatu,1), [xmin xmax],'b--','LineWidth',1.5);
fplot(@(x)fitFunctionGCompton(x,parhatu,2), [xmin xmax],'g--','LineWidth',1.5);
fplot(@(x)fitFunctionGCompton(x,parhatu,3), [xmin xmax],'c--','LineWidth',1.5);
fplot(@(x)fitFunctionGCompton(x,parhatu,4), [xmin xmax],'y--','LineWidth',1.5);
    end
else
fplot(@(x)fitFunctionLogistic(x,parhatu,1), [xmin xmax],'b--','LineWidth',1.5);
fplot(@(x)fitFunctionLogistic(x,parhatu,2), [xmin xmax],'g--','LineWidth',1.5);
fplot(@(x)fitFunctionLogistic(x,parhatu,3), [xmin xmax],'c--','LineWidth',1.5);
fplot(@(x)fitFunctionLogistic(x,parhatu,4), [xmin xmax],'y--','LineWidth',1.5);
end
b.FaceAlpha = 0.5;
xlabel('E_\gamma (11 bit ADC channels, ped sub)','FontSize',16);
ylabel(strcat('Events/channel'),'FontSize',16);
yli=get(gca,'ylim');
yl = 0.95*yli(2);
Dy = yli(2)-yli(1);
dy = 0.05*Dy;
xli=get(gca,'xlim');

for n = 1:npar
str = sprintf('%s %0.2f',parname(n),parhatu(n));
text(0.9*xli(2),yl,str);
yl = yl - dy;
hold off;
end


