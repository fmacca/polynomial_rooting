function [MSE] = mse_analytic(r,a,Sigma)
%mse_analitic computes the MSE in complex augemented notation
%from the taylor expansion
N=length(a)-1;
deriv=polyder(a);
%Jacobian matrix
delta=zeros(N,N);
for ii=1:N
    delta(ii,:)=-r(ii).^(N-1:-1:0)/polyval(deriv,r(ii));
end
MSE=[delta zeros(N,N);zeros(N,N) conj(delta)]*Sigma*[delta zeros(N,N);zeros(N,N) conj(delta)]';

end

