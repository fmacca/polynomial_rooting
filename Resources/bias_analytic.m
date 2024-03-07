function [bias] = bias_analytic(r,a,Sigma)
%bias_analytic computes the bias in complex augemented notation
%from the taylor expansion
N=length(a)-1;
deriv=polyder(a);
deriv2=polyder(deriv);
% Hessian matrices
H_aa=zeros(N,N,N); % third index refers to corresponding root
for ii=1:N
    r_=r(ii).^((N-1:-1:0)');
    r_p=((N-1:-1:0)').*r(ii).^((N-2:-1:-1)');
    H_aa(:,:,ii)=1/(polyval(deriv,r(ii)))^2*(r_*conj(r_p')+r_p*conj(r_')-polyval(deriv2,r(ii))/polyval(deriv,r(ii))*r_*conj(r_'));
end
% Computation of the bias
C=Sigma(1:N,(N+1):(2*N));
bias=zeros(2*N,1);
for ii=1:N
    bias(ii)=1/2*trace(H_aa(:,:,ii)*C);
    bias(N+ii)=conj(1/2*trace(H_aa(:,:,ii)*C));
end
end

