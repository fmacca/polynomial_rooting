clear all
close all
clc

addpath('Resources')
%% Load the dataset
% [counter r(1:N) gamma Gauss_test_HZ MBox_p Manova_d Manova_p]
N=4;
folder_name=strcat('Results/Degree0',int2str(N),'_NoiseFullMatrix_TestManova');
% folder_name=strcat('Results/Degree0',int2str(N),'_NoiseCircular_TestManova');

load(strcat(folder_name,'/','dataset.mat'));
gamma = dataset(:,2+N);
Gauss_test_HZ = dataset(:,3+N);
MBox_p = dataset(:,4+N);
Manova_d = dataset(:,5+N);
Manova_p = dataset(:,6+N);

%% Set folder for results

currDate = datestr(datetime,30);
mkdir(folder_name,'ManovaMBox');
results_folder=strcat(folder_name,'/ManovaMBox');
%All the figures we want to save should be added to the array of figures
%through figs(1)=figure(1); figs(2)=figure(2);


%% Setting the model variables
x = log(abs(gamma));
t = (MBox_p >= 0.05);

%% Logistic regression
t = categorical(t+1);
[B, ~, stats] = mnrfit(x,t); % computes the weight matrix

levels = linspace(0,1,100);
accuracies = zeros(size(levels));

for ii = 1:length(levels)
    ll = levels(ii);
    pihat = mnrval(B,x);
    t_pred = pihat(:,2)>=ll; %Predicted value
    t_pred=categorical(t_pred+1);
    cm=confusionmat(t,t_pred);

    fpr=cm(1,2)/(cm(1,2)+cm(1,1)); %1-specificity
    tpr=cm(2,2)/(cm(2,1)+cm(2,2)); %sensitivity

    accuracies(ii)=(cm(1,1)+cm(2,2))/sum(cm,"all");
end

[~,index_lev_max] = max(accuracies);
lev=levels(index_lev_max); %Model threshold

sep_line=(log((1-lev)./lev)-B(1))/B(2);
exp(sep_line)
sep_line

pihat = mnrval(B,x);
t_pred = pihat(:,2)>=lev; %Predicted value
t_pred=categorical(t_pred+1);

cm=confusionmat(t,t_pred)
fpr=cm(1,2)/(cm(1,2)+cm(1,1)); %1-specificity
tpr=cm(2,2)/(cm(2,1)+cm(2,2)); %sensitivity
accuracy=(cm(1,1)+cm(2,2))/sum(cm,"all");
fpr
tpr
accuracy

%% Plots
figs(1)=figure(1);
plot(x,MBox_p,'x'); hold on; grid on;
% plot(x,t-1,'rx');
yline(0.05,'r');
plot(x,pihat(:,2),'r.');
xline(sep_line,'b--');
legend("P-value of Box's M test","Level $\alpha=0.05$ for Box's M test","Fitted model","Separating value $\overline{\gamma}_{BoxM}$","Location","Northwest","interpreter","latex");
title("Box's M test");
xlabel("log(\gamma(z_0))");
ylabel("Probability");
hold off

figs(2)=figure(2);
plot(MBox_p,pihat(:,2),'x'); hold on; grid on;
yline(lev,'b--');
xline(0.05,'r');
xlabel("P-value of Box's M test");
ylabel("Probability from fitted model");
legend("Fitted vs actual","Chosen model threshold","Level $\alpha=0.05$ for Box's M test","Location","Southeast","interpreter","latex")
title("Logistic regression fitted probabilities vs Box's M test P-values");
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

%% Print latex table data
disp(strcat('$',num2str(sep_line),'$ & $',num2str(fpr),'$ & $',num2str(tpr),'$ & $',num2str(AUC), '$ & $',num2str(accuracy),'$'))

%% Save workspace and figures to the folder
savefig(figs,strcat(results_folder,'/figures.fig'),'compact');
clear figs
save(strcat(results_folder,'/workspace'));