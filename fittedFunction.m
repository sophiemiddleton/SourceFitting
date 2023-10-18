function [ fittedFunction ] = fittedFunction( E, par )
% fittedFunction
%   Return value of fitted spectrum at energy E
%   fcp 130333

global ME;
global NSAMPLE;
global xlo;
global xhi;

E0 = par(1);

f0=par(2);
f1=par(3);
f2 = 1-f0-f1;
A = par(4);
B = par(5);
fbg = par(6);
abg = par(7);

E1 = E0 - ME;
E2 = E1 - ME;

%A = 0.014;      %Coefficient of E(GeV)^1/4
%B = 0.03;       %Constant term


ped = 0;                    % Channels for zero energy
chPerMeV = (E0-ped)/6.13;   % Calibration constant channels/MeV
E1 = E0 - ME*chPerMeV;
E2 = E1 - ME*chPerMeV;

E0MeV = E0/chPerMeV;
E1MeV = E1/chPerMeV;
E2MeV = E2/chPerMeV;


%   Fractional resolution
s0 = sqrt(power(A/power(E0MeV/1000.,.25),2) + B^2); 
s1 = sqrt(power(A/power(E1MeV/1000.,.25),2) + B^2);
s2 = sqrt(power(A/power(E2MeV/1000.,.25),2) + B^2);
%   Absolute resolution in channels
sE0 = s0*E0MeV*chPerMeV;
sE1 = s1*E1MeV*chPerMeV;
sE2 = s2*E2MeV*chPerMeV;

norm0 = normcdf(xhi,E0,sE0) - normcdf(xlo,E0,sE0);
norm1 = normcdf(xhi,E1,sE1) - normcdf(xlo,E1,sE1);
norm2 = normcdf(xhi,E2,sE2) - normcdf(xlo,E2,sE2);
normbg = exp(-xlo/abg)-exp(-xhi/abg);



%   Fitted frequency distribution
fittedFunction = NSAMPLE*((1-fbg)*(f0*normpdf(E,E0,sE0)/norm0 + f1*normpdf(E,E1,sE1)/norm1...
    + f2*normpdf(E,E2,sE2)/norm2) + fbg*exp(-E/abg)/(abg*normbg));



end

