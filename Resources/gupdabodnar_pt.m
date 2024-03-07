function [p_value, t_value] = gupdabodnar_pt(x, Sigma0)
%NAGAO_P: This function computes the p-value 
%for Gupda and Bodnar test for the covatiance matrices
%Null hypothesis Sigma == Sigma0
%Alterantive hypothesis Sigma != Sigma0
x=x';

[p,n] = size(x);
mu_hat = sum(x,2)/n;
S = 1/(n-1) * (x-mu_hat)*(x-mu_hat)';

T = zeros(p,1);
for index=1:p
    xl0 = Sigma0(index,index);
    nu0 = Sigma0(setdiff(1:end,index),index);
    Xl0 = Sigma0(setdiff(1:end,index),setdiff(1:end,index));
    
    v = S(index,index);
    t = S(setdiff(1:end,index),index);
    W = S(setdiff(1:end,index),setdiff(1:end,index));
    
    Upsilon0 = Xl0 - nu0*nu0'/xl0;
    eta_index = sqrt(n)*Upsilon0^-0.5*(t/v-nu0/xl0)*v^0.5;
    T(index) = eta_index'*eta_index;
    
    % chi2cdf(T(index),p-1)
end

p_value = 1 - chi2cdf(max(T),p-1);
t_value = max(T);
% p_value = chi2cdf(T(1),p-1);

end