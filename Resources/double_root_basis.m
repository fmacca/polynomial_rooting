function [basis] = double_root_basis(N)
%DOUBLE_ROOT_BASIS The function returns a set of basis for double roots for
%polynomials of degree N

% We build all the possibilities af double roots
n_basis=nchoosek(N,2);
double_roots=double_root_vector_possibilities(N)';
% For each of the columns of double_roots, we complete the basis of the
% subspace
basis=zeros(N,N-1,n_basis);
for jj=1:n_basis
    basis(:,1,jj)=double_roots(:,jj);
    basis(:,2:(N-1),jj)=complete_subspace_basis(basis(:,1,jj));
end

end

