function [double_roots] = double_root_vector_possibilities(N)
%DOUBLE_ROOT_VECTOR_POSSIBILITIES Summary of this function goes here
%   Detailed explanation goes here
n_basis=nchoosek(N,2);
double_roots=zeros(n_basis,N);
jj1=1;
jj2=2;
for ii=1:n_basis
    double_roots(ii,jj1)=1;
    double_roots(ii,jj2)=1;
    if jj2 < N
        jj2=jj2+1;
    else
        jj1=jj1+1;
        jj2=jj1+1;
    end
end
    
end

