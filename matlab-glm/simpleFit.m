% simple Fit
% 
% name says it all ...
%
%

NScases = [21,28,41,51,68,73,90,110,122,127,147,170,193,207,236,262,293,310,342,373,407,428,445,474,517,547,579];
NSnewcases=[7,13,10,17,5,17,20,12,5,20,23,23,14,29,26,31,17,32,31,34,21,17,29,43,30,32];
xdata = 1:length(NScases);

ydata = [NScases,NSnewcases,NSnewcases(length(NSnewcases))];
Cinit = NScases(1);
parinit = [1.5,0.5];

parFit1 = lsqcurvefit(@(par,t) ggm(par,t,Cinit),parinit,xdata,ydata);
parFit2 = nlinfit(xdata,ydata,@(par,t) ggm(par,t,Cinit),parinit) 

ggm1 = ggm(parFit2, xdata, Cinit);
ggm1 = ggm1(1:length(NScases));

%parFit2(1)*glm1(d,rfit1,pfit1)^parFit2(2)*(1-glm1(d,rfit1,pfit1)^setalpha/setK);
subplot(2,1,1); 
scatter(xdata,NScases);
hold on;
scatter(xdata,ggm1);
subplot(2,1,2); 
scatter(1:length(NSnewcases),NSnewcases);
hold on;
%scatter(xdata,);

%ceq1:=[seq([d,glm1(d,rfit1,pfit1)],d=1..nDays)]:
%Cplot1:=plot(ceq1,style=point,color=green,symbol=solidcircle,symbolsize=16,labels=["days after March 20th","Total cases(fitted green)"],labeldirections=[horizontal,vertical]):
