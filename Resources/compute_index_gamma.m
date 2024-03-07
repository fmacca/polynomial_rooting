function [Projection_ii] = compute_index_gamma(MSE_analytic_ii, Bias_analytic_ii, r, NoiseCovarianceType)
%TODO: extend to any dimension (can be deduced from r)
    N = length(r);
    if strcmp(NoiseCovarianceType, 'circular')
        one_one = [1; 1];
        one_one_norm = abs(sqrt(one_one'*(MSE_analytic_ii(1:N,1:N)\one_one))); 
        Projection_ii= 1/one_one_norm*one_one'*(MSE_analytic_ii(1:N,1:N)\r);
    elseif strcmp(NoiseCovarianceType, 'full')
        z = [r; conj(r)];
        syms lambda
        P_inv = [[0 1; 1 0] zeros(N); zeros(N) eye(N)];
        P = inv(P_inv);
        M_lambda = (1-lambda)*MSE_analytic_ii + lambda*conj(P)'*MSE_analytic_ii*P;
        Q_lambda = (1-lambda)*eye(2*N)* + lambda * conj(P)';
        
        f = real(z'*MSE_analytic_ii*conj(Q_lambda)'*conj(M_lambda)*((MSE_analytic_ii-conj(P)'*MSE_analytic_ii*P)*conj(M_lambda)*Q_lambda-2*det(M_lambda)*(eye(2*N)-conj(P)'))*MSE_analytic_ii*z);
        solve(f)

        
    end
end