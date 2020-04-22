% GLM Fit
syms N(t) P(t) C(t);

NScases = [21,28,41,51,68,73,90,110,122,127,147,170,193,207,236,262,293,310,342,373,407,428,445,474,517,547,579];
NSnewcases=[7,13,10,17,5,17,20,12,5,20,23,23,14,29,26,31,17,32,31,34,21,17,29,43,30,32];
xdata = 1:length(NScases);

ydata = [NScases,NSnewcases,NSnewcases(length(NSnewcases))];
Cinit = NScases(1);


mu = 1/14; 
T = 300; 

parinit = [1.5,0.5,10^3,1]; % r,p,K,alpha
%parFit1 = lsqcurvefit(@(par,t) glm(par,t,Cinit),parinit,xdata,ydata);
parFit1 = nlinfit(xdata,ydata,@(par,t) glm(par,t,Cinit),parinit);
glm1 = glm(parFit1, xdata, Cinit);
glm1 = glm1(1:length(NScases));

glm1_new = (parFit1(1)*glm(parFit1, xdata, Cinit).^parFit1(2)).*((1-(glm(parFit1, xdata, Cinit)/parFit1(3)).^parFit1(4)));
glm1_new_plot = glm1_new(1:length(NSnewcases));

figure;

subplot(2,1,1); 
scatter(xdata,NScases);
hold on;
scatter(xdata,glm1);
subplot(2,1,2); 
scatter(1:length(NSnewcases),NSnewcases);
hold on;
scatter(1:length(NSnewcases),glm1_new_plot);


%Plotting N(t), I(t), and C(t) 

C_1=7;
P_1=7;
plot_spacing = linspace(1,T,300);
glm1_N_temp = (parFit1(1)*glm(parFit1, plot_spacing, Cinit).^parFit1(2)).*((1-(glm(parFit1, plot_spacing, Cinit)/parFit1(3)).^parFit1(4)));
glm1_N = glm1_N_temp(1:300);
N = spline(plot_spacing,glm1_N);
solP_temp = glm_prev_P(parFit1, plot_spacing, P_1, Cinit);
solC_temp = glm_prev_C(parFit1, plot_spacing, C_1, Cinit);

solP = solP_temp(1:300);
solC = solC_temp(1:300);

figure;
subplot(2,1,1); 
plot(plot_spacing,solC);
subplot(2,1,2); 
plot(plot_spacing, solP);
hold on;
plot(plot_spacing, ppval(N,plot_spacing));
