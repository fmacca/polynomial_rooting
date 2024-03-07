clear all
close all
clc

addpath('Resources') 
%% Parameters
N=2;%7; % order of the polynomial
SNR = [-12:10:40];%[-12:3:40];%12+[-12:3:40];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^5;%10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
K_normaltest=10^4; % Number of iterations to be used for normality test (since it cannot handle 10^5)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=10; % Number of times we generate a different polynomial

%% Generate folder for results
folder_name=strcat('Results/Degree0',int2str(N),'_NoiseFullMatrix_TestCovariance');

currDate = datestr(datetime,30);
mkdir(folder_name);
results_folder=strcat(folder_name);

rng('default')
s = rng;

%%
for counter=1:NRUNS
    % Load the dataset with previous simulations (or create a new one)
    if isfile(strcat(results_folder,'/dataset.mat'))
        load(strcat(results_folder,'/dataset.mat'));
    else
         dataset=[];
    end
    
    % Generate Polynomial and Covariance matrix
    % Generate a random circular covariance
    [Sigma,C_atilda,A] = generate_covariance(N,1,'full');
    % Generate random roots
    r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
    % Compute corresponding noise-free polynomial cefficients
    a=conj(poly(r)');

    % Simulation
    h = waitbar(0,strcat('Simulations in progress ... Please wait... ',int2str(counter),"/",int2str(NRUNS)));
    a_n=zeros(N+1,K,SNR_nsteps); %Matrix to contain coefficients at every iteration
    r_n=zeros(N,K,SNR_nsteps); %Matrix to collect roots computed at every iteration
    err_n=zeros(N,K,SNR_nsteps); %Matrix to collect the error in roots at every step
    J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
    MSE_analytic=zeros(2*N,2*N,SNR_nsteps);
    MSE_analytic_tilda=zeros(2*N,2*N,SNR_nsteps);
    Bias_analytic=zeros(2*N,SNR_nsteps);
    Bias_analytic_tilda=zeros(2*N,SNR_nsteps);
    MSE_simulated=zeros(2*N,2*N,SNR_nsteps);
    MSE_simulated_tilda=zeros(2*N,2*N,SNR_nsteps);
    Bias_simulated=zeros(2*N,SNR_nsteps);
    Bias_simulated_tilda=zeros(2*N,SNR_nsteps);
    
    % Proposed indexes of goodness
    n_basis=nchoosek(N,2);
    Projection_old=zeros(SNR_nsteps,n_basis);
    Projection_gram=zeros(SNR_nsteps,n_basis);
    Projection_from_real_projection=zeros(SNR_nsteps,n_basis);
    gamma_old=zeros(SNR_nsteps,1);
    gamma_gram=zeros(SNR_nsteps,1);
    gamma_from_real_projection=zeros(SNR_nsteps,1);

    % Normality tests
    Gauss_test_HZ=zeros(SNR_nsteps,1); % Matrices to collect the result of HZmvntest_mod

    % Tests for the mean
    HotT2_p=zeros(SNR_nsteps,1); % Hotelling T^2 test

    % Tests for the covariance
    Nagao_p=zeros(SNR_nsteps,1);
    GupdaBodnar_p=zeros(SNR_nsteps,1);
    GupdaBodnar_t=zeros(SNR_nsteps,1);

    for ii=1:SNR_nsteps
        sigma_a=(sqrt(1/SNRlin(ii)));
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
        Bias_analytic(:,ii)=bias_analytic(r,a,sigma_a^2*Sigma); % bias (complex augmented)
        Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,ii); % bias (real composite)
        
        % Projection on orthogonal of the threee spaces
        % corresponding to double roots

        % New correct version
        basis=double_root_basis(N);
        Matrix_metric_tilda=MSE_analytic_tilda(:,:,ii)+Bias_analytic_tilda(:,ii)*Bias_analytic_tilda(:,ii)';
        for bb=1:n_basis
            basis_ca=[basis(:,:,bb) 1i*basis(:,:,bb); basis(:,:,bb) -1i*basis(:,:,bb)];

            % In real coordinates
            basis_rc = 1/2*J'*basis_ca;
            orth=null((Matrix_metric_tilda\basis_rc)'); % Vector orthogonal to diagonal subspace in Malanobis metric
            orth_norm1=sqrt(orth(:,1)'*(Matrix_metric_tilda\orth(:,1))); % Norm of the orthogonal vector
            orth(:,2)=null((Matrix_metric_tilda\[basis_rc orth(:,1)])');
            orth_norm2=sqrt(orth(:,2)'*(Matrix_metric_tilda\orth(:,2)));

            Projection1=1/orth_norm1*orth(:,1)'*(Matrix_metric_tilda\([real(r); imag(r)]+Bias_analytic_tilda(:,ii)));
            Projection2=1/orth_norm2*orth(:,2)'*(Matrix_metric_tilda\([real(r); imag(r)]+Bias_analytic_tilda(:,ii)));
            tmp=(1/orth_norm1*Projection1*orth(:,1)+1/orth_norm2*Projection2*orth(:,2));
            Projection_from_real_projection(ii,bb)=tmp'*(Matrix_metric_tilda\tmp);
        end
        
        gamma_from_real_projection(ii)=min(abs(Projection_from_real_projection(ii,:)));
        
        for k=1:K
            noise_tilda=sigma_a*A*randn(2*N,1); %Generate colored noise
            a_n(:,k,ii)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
            r_curr=roots(a_n(:,k,ii)); %Compute the roots
            r_n(:,k,ii)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
            err_n(:,k,ii)=r_n(:,k,ii)-r;

            waitbar(((ii-1)*K+k)/(K*SNR_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,ii)=1/K*[err_n(:,:,ii); conj(err_n(:,:,ii))]*[err_n(:,:,ii); conj(err_n(:,:,ii))]';
        MSE_simulated_tilda(:,:,ii)=1/4*J'*MSE_simulated(:,:,ii)*J;
        Bias_simulated(:,ii)=[mean(err_n(:,:,ii),2); conj(mean(err_n(:,:,ii),2))];
        Bias_simulated_tilda(:,ii)=1/2*J'*Bias_simulated(:,ii);

        Nagao_p(ii)=nagao_p([real(r_n(:,:,ii))' imag(r_n(:,:,ii))'], MSE_analytic_tilda(:,:,ii));
        [GupdaBodnar_p(ii), GupdaBodnar_t(ii)]=gupdabodnar_pt([real(r_n(:,:,ii))' imag(r_n(:,:,ii))'], MSE_analytic_tilda(:,:,ii));
    end
    r_mean = mean(r_n,2); %Mean of the roots computed at every iteration
    
    % Save everything into a matrix [counter r1 r2 Projection Gauss_test_HZ
    % HotT2_p Nagao_p GupdaBodnar_p GupdaBodnar_t]
    dataset=[dataset; [counter*ones(SNR_nsteps,1) ones(SNR_nsteps,1)*conj(r') gamma_from_real_projection Gauss_test_HZ HotT2_p Nagao_p GupdaBodnar_p GupdaBodnar_t]];

    % Save the dataset at every iteration
    save(strcat(results_folder,'/dataset'),'dataset');
    
    close(h); %Close waitbar

end

%% Save workspace and figures to the folder
save(strcat(results_folder,'/dataset'),'dataset');