function [p_value] = nagao_p_proper_correction(x, Sigma0)
%NAGAO_P: This function computes the p-value 
%for Nagao test for the covatiance matrices
%Null hypothesis Sigma == Sigma0
%Alterantive hypothesis Sigma != Sigma0
x=x';

[p,n] = size(x);
mu_hat = sum(x,2)/n;
S = (x-mu_hat)*(x-mu_hat)';

t1 = n/2*trace(S*inv(Sigma0)/n - eye(p))^2; % sample value of Nagao's T1 statistic
% f = 1/2*p*(p+1); % degrees of freedom
f = 1/4*p*(p+2); % degrees of freedom
P_f = chi2cdf(t1,f);
P_f6 = chi2cdf(t1,f+6);
P_f4 = chi2cdf(t1,f+6);
P_f2 = chi2cdf(t1,f+6);

p_value =  1- (P_f + 1/n*(p/12*(4*p^2+9*p+7)*P_f6 - p/8*(6*p^2+13*p+9)*P_f4 + p/2*(p+1)^2*P_f2 - p/24*(2*p^2+3*p-1)*P_f));

end