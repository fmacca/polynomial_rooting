clear all
close all
clc

addpath('Resources') 
%% Generate folder for results
folder_name='Results/Degree02_NoiseCircular_VisualizingCovariance'; 

currDate = datestr(datetime,30);
mkdir(folder_name);
results_folder=strcat(folder_name);

% rng('default')
% s = rng;

%% Parameters
N=2; % order of the polynomial
SNR = [-12:3:40];
SNRlin = 10.^(SNR/10); %(sqrt(1/SNRlin(isnr)))
SNR_nsteps=numel(SNR);
K=10^4; % Number of iterations per simulation (n of noisy measurements per polynomial)
scale=2; % the roots of the polynomial will be generated with Re(r),Im(r) in [-scale +scale], a square
NRUNS=1; % Number of times we generate a different polynomial

%%
% Generate Polynomial and Covariance matrix
% Generate a random covariance
[Sigma,C_atilda,A] = generate_covariance(N,1,'circular');
% Generate random roots
r=[scale*(2*rand(N,1)-1)+scale*1i*(2*rand(N,1)-1)];
% Compute corresponding noise-free polynomial cefficients
a=conj(poly(r)');

% Simulation
h = waitbar(0,'Simulations in progress ... Please wait... ');
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

sigma_scale=zeros(SNR_nsteps,1);

% Proposed indexes of goodness
Projection1=zeros(SNR_nsteps,1);
Projection2=zeros(SNR_nsteps,1);
Projection=zeros(SNR_nsteps,1);

% Normality tests
Gauss_test_HZ=zeros(SNR_nsteps,1); % Matrices to collect the result of HZmvntest_mod

% Tests for the mean
HotT2_p=zeros(SNR_nsteps,1); % Hotelling T^2 test

for ii=1:SNR_nsteps
    sigma_a=(sqrt(1/SNRlin(ii)));
    % Compute the expected MSE matrix and bias from the analytic expression
    MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma_a^2*Sigma); % MSE matrix (complex augmented)
    MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
    Bias_analytic(:,ii)=bias_analytic(r,a,sigma_a^2*Sigma); % bias (complex augmented)
    Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,ii); % bias (real composite)
    % I do things for the projection on orthogonal of [1;1]

    sigma_scale(ii)=max(eig(sigma_a^2*Sigma));

    % Matrix_metric=MSE_analytic(:,:,ii)+Bias_analytic(:,ii)*Bias_analytic(:,ii)';
    Matrix_metric=MSE_analytic(:,:,ii);

    orth=null((Matrix_metric\[1 1i;1 1i;1 -1i;1 -1i])'); % Vector orthogonal to [1;1] in Malanobis metric
    orth_norm1=sqrt(orth(:,1)'*(Matrix_metric\orth(:,1))); % Norm of the orthogonal vector
    orth(:,2)=null((Matrix_metric\[[1 1i;1 1i;1 -1i;1 -1i] orth(:,1)])');
    orth_norm2=sqrt(orth(:,2)'*(Matrix_metric\orth(:,2)));
    
    Projection1(ii)=1/orth_norm1*orth(:,1)'*(Matrix_metric\([r; conj(r)]));
    Projection2(ii)=1/orth_norm2*orth(:,2)'*(Matrix_metric\([r; conj(r)]));
    tmp=(1/orth_norm1*Projection1(ii)*orth(:,1)+1/orth_norm2*Projection2(ii)*orth(:,2));
    Projection(ii)=tmp'*(Matrix_metric\tmp);

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
    
end
r_mean = mean(r_n,2); %Mean of the roots computed at every iteration

for ii=1:SNR_nsteps
    for k=1:K
        dist_n_from_mean(:,k,ii)=r_n(:,k,ii)-r_mean(:,:,ii);
    end
    Covariance_simulated(:,:,ii)=1/K*[dist_n_from_mean(:,:,ii); conj(dist_n_from_mean(:,:,ii))]*[dist_n_from_mean(:,:,ii); conj(dist_n_from_mean(:,:,ii))]';
    Covariance_simulated_tilda(:,:,ii)=1/4*J'*Covariance_simulated(:,:,ii)*J;
end

close(h); %Close waitbar

%% Plots
figs(1)=figure(1);

for ii=1:(2*N)
    for jj=ii:(2*N)
        if(~(ii==1 && jj==3) && ~(ii==2 && jj==4)) || ((ii==1 && jj==3) || (ii==2 && jj==4))    
            subplot(2*N,2*N,(ii-1)*2*N+jj)
            loglog(sigma_scale,abs(reshape(MSE_analytic_tilda(ii,jj,:),SNR_nsteps,1)),'-');hold on;
            loglog(sigma_scale,abs(reshape(Covariance_simulated_tilda(ii,jj,:),SNR_nsteps,1)),'xr');
            loglog(sigma_scale,abs(reshape(MSE_simulated_tilda(ii,jj,:),SNR_nsteps,1)),'xk');
            ylim([10^(-4) 10^2]);
            if(ii==1 & jj==1)
                legend([ strcat("MSE analytic"); strcat("Cov simulated"); strcat("MSE simulated")],'Location','southwest');
            end
            title(strcat("$\widetilde{MSE}_{",int2str(ii),int2str(jj),"}$ vs $\sigma$"), 'Interpreter', 'LaTeX');grid on;
        end
    end
end

%% Plots
figs(2)=figure(2);

for ii=1:(2*N)
    for jj=ii:(2*N)
        if(~(ii==1 && jj==3) && ~(ii==2 && jj==4)) || ((ii==1 && jj==3) || (ii==2 && jj==4))    
            subplot(2*N,2*N,(ii-1)*2*N+jj)
            semilogy(log(abs(Projection)),abs(reshape(MSE_analytic_tilda(ii,jj,:),SNR_nsteps,1)),'-');hold on;
            semilogy(log(abs(Projection)),abs(reshape(Covariance_simulated_tilda(ii,jj,:),SNR_nsteps,1)),'xr');
            semilogy(log(abs(Projection)),abs(reshape(MSE_simulated_tilda(ii,jj,:),SNR_nsteps,1)),'xk');
            ylim([10^(-4) 10^2]);
            if(ii==1 & jj==1)
                legend([ strcat("MSE analytic"); strcat("Cov simulated"); strcat("MSE simulated")],'Location','southwest');
            end
            title(strcat("$\widetilde{MSE}_{",int2str(ii),int2str(jj),"}$ vs $\sigma$"), 'Interpreter', 'LaTeX');grid on;
        end
    end
end

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));