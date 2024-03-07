function ord_index = order_roots_permutations_hermitian(r,centers)
%Order_roots: receives N roots and N centers,
%reorders the roots by selecting the permutation 
%that mimimizes the distance.
N = length(r);
P = perms(1:N);
n_perm=size(P,1);
dist=zeros(n_perm,1);

for ii=1:n_perm
    dist(ii)=(r(P(ii,:))-centers)'*(r(P(ii,:))-centers);
end
[~,indx]=min(dist);
ord_index = P(indx,:);

end

