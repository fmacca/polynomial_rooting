function [Sigma,C_atilda,A] = generate_covariance(N,sigma_a,type)
%generate_covariance
%this is a function to support the simulations.
%It returns the matrix complex augmented notation, in real composite notation and its cholesky decomposition
%type can be either 'full' or 'circular'

flag=0; %flag to become 1 when a valid (spd) matrix is generated
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)];
while(~flag)
    Temp = sigma_a*1/2/N*wishrnd(eye(2*N),2*N);
    Sigma=J*Temp*J';
    if(strcmp(type,'circular'))
        Gamma=Sigma(1:N,1:N);
        C_atilda=1/4*J'*[Gamma zeros(N);zeros(N) conj(Gamma)]*J;
    elseif(strcmp(type,'full'))
        C_atilda=1/4*J'*Sigma*J;
    else
        disp('Unsupported type parameter inserted')
    end
    try A=chol(C_atilda)';
        flag=1;
    catch ME
        disp('Attempted matrix is not symmetric positive definite')
    end
    Sigma=J*C_atilda*J';   
end

end

