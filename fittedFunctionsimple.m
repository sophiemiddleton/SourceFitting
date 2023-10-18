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
sigma = par(4);
fbg = par(5);
abg = par(6);


E1 = E0 - ME;
E2 = E1 - ME;

%A = 0.014;      %Coefficient of E(GeV)^1/4
%B = 0.03;       %Constant term

%   Fractional resolution
s0 = sqrt(power(A/power(E0/1000.,.25),2) + B^2); 
s1 = sqrt(power(A/power(E1/1000.,.25),2) + B^2);
s2 = sqrt(power(A/power(E2/1000.,.25),2) + B^2);
%   Absolute resolution
sE0 = s0*E0;
sE1 = s1*E1;
sE2 = s2*E2;


norm0 = normcdf(xhi,E0,sE0) - normcdf(xlo,E0,sE0);
norm1 = normcdf(xhi,E1,sE1) - normcdf(xlo,E1,sE1);
norm2 = normcdf(xhi,E2,sE2) - normcdf(xlo,E2,sE2);
normbg = exp(-xlo/abg)-exp(-xhi/abg);



%   Fitted frequency distribution
fittedFunction = NSAMPLE*((1-fbg)*(f0*normpdf(E,E0,sE0)/norm0 + f1*normpdf(E,E1,sE1)/norm1...
    + f2*normpdf(E,E2,sE2)/norm2) + fbg*exp(-E/abg)/(abg*normbg));



end

