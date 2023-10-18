% Normal moment calculator
% Untested, unfinished
% Reference: https://en.wikipedia.org/wiki/Normal_distribution
% fcp 170920 
p = -1;
a = -p/2;
b = 1/2;
z = 1;
Mabz = hypergeom(a,b,z);
M2abz = hypergeom(a+1-b,2-b,z);
A = Mabz*gamma(1-b)/gamma(a+1-b);
A2 = M2abz*gamma(b-1)/gamma(a);
A2z = A2*z^(1-b);
