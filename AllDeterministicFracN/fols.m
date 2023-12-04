function [b, stderr, tstat, Var_b, s2, yhat, Rsquare] = fols(X,y)
% OLS regression
n = size(y,1); k = size(X,2);
XtX = X'*X; 
invXtX = XtX \ eye(k); 
b = invXtX *X'*y;
yhat = X*b;
e = y-yhat; 
RSS = e'*e;
s2 = RSS/(n-k);
Var_b = s2 * invXtX;
stderr = sqrt(diag(Var_b));
tstat = b ./ stderr;
TSS = sum( (y - mean(y)).^2);
Rsquare = 1- RSS/TSS;

end
