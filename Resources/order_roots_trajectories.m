function ord_index = order_roots_trajectories(centers,a,a_mod)
%Order_roots: receives N roots and N centers,
%reorders the roots by following a trajectory.
%At each step it selects the permutation 
%that mimimizes the distance.

scale_factor=5;
%compute minimum distance between roots
n_steps=10;%should be based on distances

lambdas=linspace(0,1,n_steps);
for kk=2:n_steps
    lambda=lambdas(kk);
    a_curr=lambda*a_mod+(1-lambda)*a;
    r_curr=roots(a_curr);
    ord_index=order_roots_permutations(r_curr,centers);
    centers=r_curr(ord_index); % Order the estimated roots
end

end

