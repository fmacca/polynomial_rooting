function max_abs_dist = max_abs_roots_sets(r,centers)
%TODO
N = length(r);
P = perms(1:N);
n_perm=size(P,1);
dist=zeros(n_perm,1);

for ii=1:n_perm
    dist(ii)=max(abs(r(P(ii,:))-centers));
end
max_abs_dist=min(dist);

end

