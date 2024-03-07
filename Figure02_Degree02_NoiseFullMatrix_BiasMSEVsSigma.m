%% Bias vs Sigma Figure 1 Special case, but with Non-Proper noise
clear all
close all
clc

addpath('Resources') 

%% Generate folder for results
folder_name='Results/Degree02_NoiseFullMatrix_BiasMSEVsSigma'; %Name for the results folder: it should be named after the kind of test performed

currDate = datestr(datetime,30);
mkdir(folder_name,currDate);
results_folder=strcat(folder_name,'/',currDate);

%% Parameters
N=2; % order of the polynomial
sigma_scale = logspace(-3,0,40);%[0.01:0.01:0.4]; % standard deviation of the noise on coefficients
sigma_nsteps=numel(sigma_scale); 
K=10^4;%10^6; % Number of iterations per simulation (n of noisy measurements per polynomial)
Repeat_experiment = 2; % 10^2

%% Generate Polynomial and Covariance matrix
J=[eye(N) 1i*eye(N);eye(N) -1i*eye(N)]; % Notation change matrix
% Generate a random circular covariance
%Sigma = toeplitz([1 0.5 0.25 0.125]); % We randomize here
Sigma = [eye(N) 0.5*1i*eye(N);-0.5*1i*eye(N) eye(N)];
%Sigma = Sigma/max(eig(Sigma)); % We normalize for comparison!
C_atilda=1/4*J'*Sigma*J;
A=chol(C_atilda)';
% Generate random roots
r=[1; 1i];
% Compute corresponding noise-free polynomial cefficients
a=conj(poly(r)');

%% Simulation
a_n=zeros(N+1,K,sigma_nsteps,Repeat_experiment); %Matrix to contain coefficients at every iteration
r_n=zeros(N,K,sigma_nsteps,Repeat_experiment); %Matrix to collect roots computed at every iteration
err_n=zeros(N,K,sigma_nsteps,Repeat_experiment); %Matrix to collect the error in roots at every step

MSE_analytic=zeros(2*N,2*N,sigma_nsteps); % Repeat_experiment dimension not added here, because this does not change!
MSE_analytic_tilda=zeros(2*N,2*N,sigma_nsteps);
Bias_analytic=zeros(2*N,sigma_nsteps); % Repeat_experiment dimension not added here, because this does not change!
Bias_analytic_tilda=zeros(2*N,sigma_nsteps);
MSE_simulated=zeros(2*N,2*N,sigma_nsteps,Repeat_experiment); % The MSE of the sample
MSE_simulated_tilda=zeros(2*N,2*N,sigma_nsteps,Repeat_experiment);
Bias_simulated=zeros(2*N,sigma_nsteps,Repeat_experiment); % The MSE of the sample
Bias_simulated_tilda=zeros(2*N,sigma_nsteps,Repeat_experiment);
MSE_estimator_simulated=zeros(2*N,2*N,sigma_nsteps); % The MSE of the estimator for the roots (ie the average root)
MSE_estimator_simulated_tilda=zeros(2*N,2*N,sigma_nsteps);
for rep=1:Repeat_experiment % TODO: convert this to something that saves and loads averages to disk
    h = waitbar(0,strcat('Simulation ',int2str(rep),'/',int2str(Repeat_experiment),' in progress ... Please wait...'));
    for ii=1:sigma_nsteps
        sigma=sigma_scale(ii);
        % Compute the expected MSE matrix and bias from the analytic expression
        MSE_analytic(:,:,ii)=mse_analytic(r,a,sigma^2*Sigma); % MSE matrix (complex augmented)
        MSE_analytic_tilda(:,:,ii)=1/4*J'*MSE_analytic(:,:,ii)*J; % MSE matrix (real composite)
        Bias_analytic(:,ii)=bias_analytic(r,a,sigma^2*Sigma); % bias (complex augmented)
        Bias_analytic_tilda(:,ii)=1/2*J'*Bias_analytic(:,ii); % bias (real composite)
        
        for k=1:K
            noise_tilda=sigma*A*randn(2*N,1); %Generate colored noise
            a_n(:,k,ii,rep)=a+[0;noise_tilda(1:N)+1i*noise_tilda(N+1:2*N)]; %Add noise to coefficients
            r_curr=roots(a_n(:,k,ii,rep)); %Compute the roots
            r_n(:,k,ii,rep)=r_curr(order_roots_permutations(r_curr,r)); %Save roots ordered w.r.t. original roots
            err_n(:,k,ii,rep)=r_n(:,k,ii,rep)-r;
    
            waitbar(((ii-1)*K+k)/(K*sigma_nsteps)) %Update waitbar
        end
        MSE_simulated(:,:,ii,rep)=1/K*[err_n(:,:,ii,rep); conj(err_n(:,:,ii,rep))]*[err_n(:,:,ii,rep); conj(err_n(:,:,ii,rep))]';
        MSE_simulated_tilda(:,:,ii,rep)=1/4*J'*MSE_simulated(:,:,ii,rep)*J;
        Bias_simulated(:,ii,rep)=[mean(err_n(:,:,ii,rep),2); conj(mean(err_n(:,:,ii,2),2))];
        Bias_simulated_tilda(:,ii,rep)=J'*Bias_simulated(:,ii,rep); % TODO: IT SHOULD NOT BE COMPLEX
    end
    
    close(h); %Close waitbar
end

% Compute the MSE of the estimator
r_mean = reshape(mean(r_n,2),N,sigma_nsteps,Repeat_experiment); %Mean of the roots computed at every iteration
err_mean = zeros(N,sigma_nsteps,Repeat_experiment);
for ii = 1:sigma_nsteps
    for rep = 1:Repeat_experiment
        err_mean(:,ii,rep) = r_mean(:,ii,rep) - r;
    end
    err_mean_ii = reshape(err_mean(:,ii,:),N,Repeat_experiment);
    MSE_estimator_simulated(:,:,ii) = 1/Repeat_experiment*[err_mean_ii; conj(err_mean_ii)]*[err_mean_ii; conj(err_mean_ii)]';
    MSE_estimator_simulated_tilda(:,:,ii) = 1/4*J'*MSE_estimator_simulated(:,:,ii)*J;
end


%% Theoretical bound for the bias
% (pimi 28 valori sono accettabili, ma i primi sono poco affidabili)
% (limite su sigma per validita del bound 0.280262)
theoretical_boud_sigmas = [0.01:0.01:0.4];
theoretical_bound_bias = [1.05833*10^-12, 2.64577*10^-13, 1.17586*10^-13, ...
    6.61346*10^-14, 4.23273*10^-14, 2.9407*10^-14, 2.08428*10^-12, ...
    2.92368*10^-9, 4.03287*10^-7, 0.0000131708, 0.000168481, 0.00114202, ...
    0.00496111, 0.0156428, 0.0389431, 0.0811598, 0.14762, 0.241553, 0.363658, ...
    0.512309, 0.68413, 0.874691, 1.07914, 1.2927, 1.51099, 1.73022, ...
    1.94724, 2.15958, 2.36534, 2.56319, 2.75222, 2.93189, 3.10197, ...
    3.26245, 3.41347, 3.55531, 3.68836, 3.81301, 3.92974, 4.039];
theoretical_bound_sigmas_to_plot = theoretical_boud_sigmas(1:28);
theoretical_bound_bias_to_plot = theoretical_bound_bias(1:28);

%% We want to plot the value of gamma
% Proposed indexes of goodness
Projection1=zeros(sigma_nsteps,1);
Projection2=zeros(sigma_nsteps,1);
Projection=zeros(sigma_nsteps,1);

for ii=1:sigma_nsteps
    Matrix_metric=MSE_analytic(:,:,ii);
    
    orth=null((Matrix_metric\[1 1i;1 1i;1 -1i;1 -1i])'); % Vector orthogonal to [1;1] in Malanobis metric
    orth_norm1=sqrt(orth(:,1)'*(Matrix_metric\orth(:,1))); % Norm of the orthogonal vector
    orth(:,2)=null((Matrix_metric\[[1 1i;1 1i;1 -1i;1 -1i] orth(:,1)])');
    orth_norm2=sqrt(orth(:,2)'*(Matrix_metric\orth(:,2)));
    
    % Projection1(ii)=1/orth_norm1*orth(:,1)'*(Matrix_metric\([r; conj(r)]+Bias_analytic(:,ii)));
    % Projection2(ii)=1/orth_norm2*orth(:,2)'*(Matrix_metric\([r; conj(r)]+Bias_analytic(:,ii)));
    Projection1(ii)=1/orth_norm1*orth(:,1)'*(Matrix_metric\([r; conj(r)]));
    Projection2(ii)=1/orth_norm2*orth(:,2)'*(Matrix_metric\([r; conj(r)]));
    tmp=(1/orth_norm1*Projection1(ii)*orth(:,1)+1/orth_norm2*Projection2(ii)*orth(:,2));
    Projection(ii)=tmp'*(Matrix_metric\tmp);
end


%% Estimator MSE Complex Representation

figs(1)=figure(1);
subplot(1,1,1)
loglog(theoretical_bound_sigmas_to_plot, theoretical_bound_bias_to_plot.^2, 'r--');hold on;grid on;

loglog(sigma_scale, reshape((abs(MSE_analytic(1,1,:))+abs(MSE_analytic(2,2,:)))/K,sigma_nsteps,1), 'bx');hold on;

loglog(sigma_scale, reshape((abs(MSE_estimator_simulated(1,1,:))+abs(MSE_estimator_simulated(2,2,:))),sigma_nsteps,1),'k-','LineWidth',2);hold on;

axis([sigma_scale(1) sigma_scale(end) 10^-9 10^-1])

title("MSE of the estimator vs \sigma")

legend(["Theoretical threshold bias", "Analytical MSE prediction", "MSE from simualation"],'Location','northwest')


%% Estimator MSE Real Representation

figs(2)=figure(2);
subplot(1,1,1)

% loglog(theoretical_bound_sigmas_to_plot, theoretical_bound_bias_to_plot.^2, 'r--');hold on;

toplot = abs(MSE_analytic_tilda(1,1,:));
for iiii=2:4
    toplot = toplot + abs(MSE_analytic_tilda(iiii,iiii,:));
end
analytical_MSE = reshape(toplot,sigma_nsteps,1)/K;
loglog(sigma_scale, analytical_MSE, 'b-');hold on;grid on;

toplot = abs(MSE_estimator_simulated_tilda(1,1,:));
for iiii=2:4
    toplot = toplot + abs(MSE_estimator_simulated_tilda(iiii,iiii,:));
end
loglog(sigma_scale, reshape(toplot,sigma_nsteps,1), 'kx','LineWidth',2);hold on;

% toplot = abs(sum(Bias_analytic_tilda(1,:,:),3));
% for iiii=2:4
%     toplot = toplot + abs(sum(Bias_analytic_tilda(iiii,:,:),3));
% end
% loglog(sigma_scale, abs(toplot.^2), 'b--','LineWidth',2);hold on;


axis([sigma_scale(1) sigma_scale(end) 10^-11 10^-1])

title("MSE of the estimator vs \sigma")

legend(["Analytical MSE prediction", "MSE from simualation"],'Location','northwest')

%% Two clouds of points
figs(3) = figure(3);

d = 10; % A number in 1:sigma_nsteps
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n(ii,:,d,rep)),imag(r_n(ii,:,d,rep)),'.','MarkerSize',1); hold on; % Simulated roots
end
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii],d)),[real(r(ii)),imag(r(ii))])
% end
plot(real(r_mean(:,d,rep)),imag(r_mean(:,d,rep)),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(2*[-1,1,-1,1]);
title(strcat("Distribution of the roots with \sigma = ",num2str(sigma_scale(d))));grid on;hold off

figs(4) = figure(4);

d = 39; % A number in 1:sigma_nsteps
viscircles([0 0],1,'color','b','linestyle','--','LineWidth',0.1);hold on;
plot(zeros(2,1),5*[-1,1],'b--','LineWidth',0.1);plot(5*[-1,1],zeros(2,1),'b--','LineWidth',0.1);
for ii=1:N
    plot(real(r_n(ii,:,d,rep)),imag(r_n(ii,:,d,rep)),'.','MarkerSize',1); hold on; % Simulated roots
end
% for ii=1:N
% %     Ellipse_plot(0.1*inv(Cov_ztilda([ii N+ii],[ii N+ii])),[real(avgs(ii))-real(bias(ii)),imag(avgs(ii))-imag(bias(ii))])
%     ellipse_plot(0.1*inv(MSE_analytic_tilda([ii N+ii],[ii N+ii],d)),[real(r(ii)),imag(r(ii))])
% end
plot(real(r_mean(:,d,rep)),imag(r_mean(:,d,rep)),'.b','MarkerSize',15); % Mean of estimated roots
plot(real(r),imag(r),'*k','MarkerSize',20); % True roots
axis equal;axis(2*[-1,1,-1,1]);
title(strcat("Distribution of the roots with \sigma =",num2str(sigma_scale(d))));grid on;hold off

%% Estimator Bias Real Representation

figs(5)=figure(5);
subplot(1,1,1)


toplot = abs(sum(Bias_analytic_tilda(1,:,:),3));
for iiii=2:4
    toplot = toplot + abs(sum(Bias_analytic_tilda(iiii,:,:),3));
end
loglog(sigma_scale, abs(toplot), 'b-');hold on;grid on;

% Compute the bias of the estimator
r_mean = reshape(mean(r_n,2),N,sigma_nsteps,Repeat_experiment); %Mean of the roots computed at every iteration
err_mean = zeros(N,sigma_nsteps,Repeat_experiment);
Bias_estimator_simulated=zeros(2*N,sigma_nsteps); % The bias of the estimator for the roots (ie the average root)
Bias_estimator_simulated_tilda=zeros(2*N,sigma_nsteps);
for ii = 1:sigma_nsteps
    for rep = 1:Repeat_experiment
        err_mean(:,ii,rep) = r_mean(:,ii,rep) - r;
    end
    err_mean_ii = reshape(err_mean(:,ii,:),N,Repeat_experiment);
    Bias_estimator_simulated(:,ii) = 1/Repeat_experiment*sum([err_mean_ii; conj(err_mean_ii)],2);
    Bias_estimator_simulated_tilda(:,ii) = 1/2*J'*Bias_estimator_simulated(:,ii);
end


toplot = abs(Bias_estimator_simulated_tilda(1,:));
for iiii=2:4
    toplot = toplot + abs(Bias_estimator_simulated_tilda(iiii,:));
end
loglog(sigma_scale, reshape(toplot,sigma_nsteps,1), 'kx','LineWidth',2);hold on;


axis([sigma_scale(1) sigma_scale(end) 10^-7 10^0])

title("Bias of the estimator vs \sigma")

legend(["Analytical Bias prediction", "Bias from simualation"],'Location','northwest')

%% Gamma

fig(6) = figure(6);

subplot(5,1,[1 2])
% loglog(theoretical_bound_sigmas_to_plot, theoretical_bound_bias_to_plot.^2, 'r--');hold on;

toplot = abs(MSE_analytic_tilda(1,1,:));
for iiii=2:4
    toplot = toplot + abs(MSE_analytic_tilda(iiii,iiii,:));
end
analytical_MSE = reshape(toplot,sigma_nsteps,1)/K;
loglog(sigma_scale, analytical_MSE, 'b-');hold on;grid on;

toplot = abs(MSE_estimator_simulated_tilda(1,1,:));
for iiii=2:4
    toplot = toplot + abs(MSE_estimator_simulated_tilda(iiii,iiii,:));
end
loglog(sigma_scale, reshape(toplot,sigma_nsteps,1), 'kx','LineWidth',2);hold on;

% toplot = abs(sum(Bias_analytic_tilda(1,:,:),3));
% for iiii=2:4
%     toplot = toplot + abs(sum(Bias_analytic_tilda(iiii,:,:),3));
% end
% loglog(sigma_scale, abs(toplot.^2), 'b--','LineWidth',2);hold on;

axis([sigma_scale(1) sigma_scale(end) 10^-11 10^-1])
title("MSE of the estimator vs \sigma")
legend(["Analytical MSE prediction", "MSE from simualation"],'Location','northwest')

subplot(5,1,[3 4])
toplot = abs(sum(Bias_analytic_tilda(1,:,:),3));
for iiii=2:4
    toplot = toplot + abs(sum(Bias_analytic_tilda(iiii,:,:),3));
end
loglog(sigma_scale, abs(toplot), 'b-');hold on;grid on;

% Compute the bias of the estimator
r_mean = reshape(mean(r_n,2),N,sigma_nsteps,Repeat_experiment); %Mean of the roots computed at every iteration
err_mean = zeros(N,sigma_nsteps,Repeat_experiment);
Bias_estimator_simulated=zeros(2*N,sigma_nsteps); % The bias of the estimator for the roots (ie the average root)
Bias_estimator_simulated_tilda=zeros(2*N,sigma_nsteps);
for ii = 1:sigma_nsteps
    for rep = 1:Repeat_experiment
        err_mean(:,ii,rep) = r_mean(:,ii,rep) - r;
    end
    err_mean_ii = reshape(err_mean(:,ii,:),N,Repeat_experiment);
    Bias_estimator_simulated(:,ii) = 1/Repeat_experiment*sum([err_mean_ii; conj(err_mean_ii)],2);
    Bias_estimator_simulated_tilda(:,ii) = 1/2*J'*Bias_estimator_simulated(:,ii);
end

toplot = abs(Bias_estimator_simulated_tilda(1,:));
for iiii=2:4
    toplot = toplot + abs(Bias_estimator_simulated_tilda(iiii,:));
end
loglog(sigma_scale, reshape(toplot,sigma_nsteps,1), 'kx','LineWidth',2);hold on;

axis([sigma_scale(1) sigma_scale(end) 10^-7 10^0])
title("Bias of the estimator vs \sigma")
legend(["Analytical Bias prediction", "Bias from simualation"],'Location','northwest')

subplot(5,1,5)
semilogx(sigma_scale,log(abs(Projection))); hold on; grid on;
threshold_value = 5.10;
yline(threshold_value,'r--');
title("Gamma at level \sigma vs threshold for \gamma_{T^2}")
legend('\gamma computed at \sigma','Threshold \gamma_{T^2}')


%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));
