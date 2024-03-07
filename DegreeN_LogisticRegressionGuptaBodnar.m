clear all
close all
clc

addpath('Resources')
%% Load the dataset
% [counter r1 r2 Projection Gauss_test_HZ HotT2_p Nagao_p GupdaBodnar_p GupdaBodnar_t]
N=3;
folder_name=strcat('Results/Degree0',int2str(N),'_NoiseFullMatrix_TestCovariance');
% folder_name=strcat('Results/Degree0',int2str(N),'_NoiseCircular_TestCovariance');

load(strcat(folder_name,'/','dataset.mat'));
gamma = dataset(:,2+N);
Gauss_test_HZ = dataset(:,3+N);
HotT2_p = dataset(:,4+N);
Nagao_p = dataset(:,5+N);
GupdaBodnar_p = dataset(:,6+N);
GupdaBodnar_t = dataset(:,7+N);

%% Set folder for results

currDate = datestr(datetime,30);
mkdir(folder_name,'GupdaBodnar');
results_folder=strcat(folder_name,'/GupdaBodnar');
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);


%% Setting the model variables
x = log(abs(gamma));
t = (GupdaBodnar_t <= 10.25);

%% Logistic regression
t = categorical(t+1);
[B, ~, stats] = mnrfit(x,t); % computes the weight matrix

% Model threshold
lev=0.4;
sep_line=(log((1-lev)./lev)-B(1))/B(2);
exp(sep_line)
sep_line

pihat = mnrval(B,x);
t_pred = pihat(:,2)>=lev; %Predicted value
t_pred=categorical(t_pred+1);

cm=confusionmat(t,t_pred)
fpr=cm(1,2)/(cm(1,2)+cm(1,1)); %1-specificity
tpr=cm(2,2)/(cm(2,1)+cm(2,2)); %sensitivity
fpr
tpr


%% Plots
figs(1)=figure(1);
plot(x,GupdaBodnar_p,'x'); hold on; grid on;
% plot(x,t-1,'rx');
yline(0.05,'r');
plot(x,pihat(:,2),'r.');
xline(sep_line,'b--');
legend("P-value of GupdaBodnar test","Level $\alpha=0.05$ for GupdaBodnar test","Fitted model","Separating value $\overline{\gamma}_{GupdaBodnar}$","Location","Northwest","interpreter","latex");
title("GupdaBodnar test");
xlabel("log(|\gamma(z_0)|)");
ylabel("Probability");
hold off

figs(2)=figure(2);
plot(GupdaBodnar_p,pihat(:,2),'x'); hold on; grid on;
yline(lev,'b--');
xline(0.05,'r');
xlabel("P-value of GupdaBodnar test");
ylabel("Probability from fitted model");
legend("Fitted vs actual","Chosen model threshold","Level $\alpha=0.05$ for GupdaBodnar test","Location","Southeast","interpreter","latex")
title("Logistic regression fitted probabilities vs GupdaBodnar test P-values");
hold off

figs(3)=figure(3); % ROC curve
% [X,Y]=perfcurve(t,x,2);
[X,Y,T,AUC] = perfcurve(t,x,2);
AUC
plot(X,Y);
xlabel('False positive rate') 
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression')
hold on; grid on
plot(fpr,tpr,'rx');

figs(4)=figure(4);
semilogy(x,GupdaBodnar_t,'x'); hold on; grid on;
% plot(x,t-1,'rx');
yline(gupdabodnar_thesholds(4,0.05),'r');
% plot(x,pihat(:,2),'r.');
xline(sep_line,'b--');
legend("T-value of GupdaBodnar test","Threshold","Separating value $\overline{\gamma}_{GupdaBodnar}$","Location","Northwest","interpreter","latex");
title("GupdaBodnar test");
xlabel("log(|\gamma(z_0)|)");
ylabel("T value");
hold off

disp(strcat('$',num2str(sep_line),'$ & $',num2str(fpr),'$ & $',num2str(tpr),'$ & $',num2str(AUC),'$'))


%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));