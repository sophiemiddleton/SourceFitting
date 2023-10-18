function [ m2lnL ] = likelihoodE0( par )
% likelihoodE0
%   Return -2 ln L for energy spectrum at location parameter E0
%   fcp 170617

global ME;
global NSAMPLE;
global EDATA;
global xlo;     % low end of fit range (channels)
global xhi;     % high end of fit range (channels)

E0 = par(1);    % Energy of primary photon, in channels
f0 = par(2);    % Fraction of signal in primary photon
f1 = par(3);    % Fraction of signal in first escape peak
A = par(4);     % E^1/4 coefficient
B = par(5);     % Constant term in Eres 
fbg = par(6);   % Fraction in background
abg = par(7);   % Exponential slope

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

f2 = 1-f0-f1;

Abg = 1/(abg*(exp(-xlo/abg) - exp(-xhi/abg)));
Asig = 1/sqrt(2*pi);

m2lnL = 0;


for i = 1:NSAMPLE;
    E = EDATA(i);
    m2lnL = m2lnL + log( ...
        (1-fbg)*Asig*( ...
            (f0/sE0)*exp(-.5*((E-E0)/sE0)^2)...
            + (f1/sE1)*exp(-.5*((E-E1)/sE1)^2) ...
            + (f2/sE2)*exp(-.5*((E-E2)/sE2)^2) ...
            )...
        + fbg*exp(-E/abg)*Abg ...
        );
end

m2lnL = -2.*m2lnL;

end

