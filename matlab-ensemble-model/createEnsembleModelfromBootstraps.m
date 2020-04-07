clear 
close all

addpath('c:/Users/colem/Documents/GitHub/mcmcstat')

% Code to generate ensemble forecast based on the paper:
% Real-time forecasting of epidemic trajectories using computational dynamic ensembles
% G. Chowell1,2,*, R. Luo1, K. Sun2, K. Roosa1, A. Tariq1, C. Viboud2

%RMSEC contains the RMSE of the calibration period for model 1 (GGM) and
%model 2 (GLM)

RMSEC=[15.0582   18.7016];

weights=1./RMSEC;

%Normalize weights
weights=weights./sum(weights);

%load model datasets to construct the ensemble model
load 'models-data.mat'

%forecastss1 contains M mean fits for model 1 (M columns)
%forecastss2 contains M mean fits for model 2 (M columns)

% Lines 25-39 generate the ensemble model curves (A total of MxM curves from our two component models)
forecastss=[]; % in this variable we save the new ensemble model curves

for t=1:length(forecastss1(:,1))
    t
    
    forecastss3=[];
    
    for i=1:M
        forecastss3=[forecastss3 [weights(1)*forecastss1(t,i)+weights(2)*forecastss2(t,:)]];
        
    end
    
    forecastss=[forecastss;forecastss3];
    
end

%ensemble model is now generated

figure

timevect=(1:length(data(:,1)))'*DT;

% ENSEMBLE forecast

%subplot(1,3,3)


line1=plot(timevect2,plims(forecastss',0.5),'r');
set(line1,'LineWidth',2)
hold on



line1=plot(timevect2,plims(forecastss',0.025),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect2,plims(forecastss',0.975),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect,data(:,2),'bo')
set(line1,'LineWidth',2,'Markersize',6)

xlabel('\fontsize{24}Time (days)');
ylabel('\fontsize{24}Case incidence')

axis([0 (50+forecastingperiod)*DT 0 max(data(:,2)+50)])
%axis([timevect2(1) data(tf2+forecastingperiod+10,1)*DT 0 max(data(:,2))*2])

line2=[t_window(end)*DT 0;t_window(end)*DT max(data(:,2)+200)];

%line2=[timevect(end) 0;timevect(end) max(data(:,2)+200)];

plot(line2(:,1),line2(:,2),'k--')

set(gca,'FontSize', 16);
set(gcf,'color','white')



% Model 1 (GGM) forecast
%subplot(1,3,1)
figure

line1=plot(timevect2,plims(forecastss1',0.5),'r')
set(line1,'LineWidth',2)
hold on

line1=plot(timevect2,plims(forecastss1',0.025),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect2,plims(forecastss1',0.975),'r--')
set(line1,'LineWidth',2)


line1=plot(timevect,data(:,2),'bo')
set(line1,'LineWidth',2,'Markersize',6)

xlabel('\fontsize{24}Time (days)');
ylabel('\fontsize{24}Case incidence')

axis([0 (50+forecastingperiod)*DT 0 max(data(:,2)+50)])

line2=[t_window(end)*DT 0;t_window(end)*DT max(data(:,2)+200)];

plot(line2(:,1),line2(:,2),'k--')

set(gca,'FontSize', 16);
set(gcf,'color','white')




% Model 2 (GRM forecast)

%subplot(1,3,2)
figure

line1=plot(timevect2,plims(forecastss2',0.5),'r')
set(line1,'LineWidth',2)
hold on

line1=plot(timevect2,plims(forecastss2',0.025),'r--')
set(line1,'LineWidth',2)

line1=plot(timevect2,plims(forecastss2',0.975),'r--')
set(line1,'LineWidth',2)


line1=plot(timevect,data(:,2),'bo')
set(line1,'LineWidth',2,'Markersize',6)

xlabel('\fontsize{24}Time (days)');
ylabel('\fontsize{24}Case incidence')

axis([0 (50+forecastingperiod)*DT 0 max(data(:,2)+50)])
%axis([timevect2(1) data(tf2+forecastingperiod+10,1)*DT 0 max(data(:,2))*2])

line2=[t_window(end)*DT 0;t_window(end)*DT max(data(:,2)+200)];

%line2=[timevect(end) 0;timevect(end) max(data(:,2)+200)];

plot(line2(:,1),line2(:,2),'k--')

set(gca,'FontSize', 16);
set(gcf,'color','white')


