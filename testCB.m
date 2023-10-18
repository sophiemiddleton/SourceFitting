% Script testCB
% Test CB functions
% fcp 170629

TESTfunctionCB = true;
TESTfunctionCBintegral = true;

if(TESTfunctionCB)
    n = 3.;
    mu = 0.;
    sigma = 1.;
    alpha = 1.5;
    xinterval = [-5,5];
    f = @(x)functionCB(x,alpha,n,mu,sigma);
    fplot(f, xinterval)
    x = -1:.05:5;
    xt = x';
    y=functionCB(xt,alpha,n,mu,sigma);
end
   
if(TESTfunctionCBintegral)
    n = 2.;
    mu = 0.;
    sigma = 1.;
    alpha = 1.;
    xinterval = [-5,5];
    f = @(x)functionCBintegral(x,alpha,n,mu,sigma);
    fplot(f, xinterval)
    hold on;
    alpha = -alpha;
    f = @(x)functionCBintegral(x,alpha,n,mu,sigma);
    fplot(f, xinterval)
    hold off;
    x = -1:.05:5;
    xt = x';
    y=functionCB(xt,alpha,n,mu,sigma);
end