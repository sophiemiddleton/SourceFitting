% Script testLogistic
% Test Logistic functions
% fcp 170630

TESTfunctionLogistic = true;
TESTfunctionLogisticintegral = true;

if(TESTfunctionLogistic)
    beta = 0.5;
    x0 = 1.;
    xinterval = [-5,5];
    f = @(x)functionLogistic(x,x0,beta);
    figure;
    fplot(f, xinterval)
    title('TESTfunctionLogistic');
%    x = -1:.05:5;
%    xt = x';
%    y=functionLogistic(xt,beta);
end
   
if(TESTfunctionLogisticintegral)
    beta = 0.5;
    x0 = 1.;
    xlow = -2.;
    xinterval = [xlow,5];
    f = @(x)functionLogisticintegral(xlow,x,x0,beta);
    figure;
    fplot(f, xinterval)
    title('TESTfunctionLogisticintegral');
end