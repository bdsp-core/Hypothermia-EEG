clear all; clc; format compact; 


td = '/Users/mbw/cdac Dropbox/Datasets/zz_Hypothermia/Hypothermia-EEG-Data/'; 
full_path = [td 'BSPTEMP']
% load BSPTEMP % from FIG8A
n=1; 
x=x(1:n:end); 
y=y(1:n:end); 

yy=log(y./(1-y));
[b,dev,stats]=glmfit(x,yy,'normal');
c=corr(x',yy')
ci = bootci(5000,@corr,x',yy')
xx=linspace(min(x),max(x),100); 
yfit = glmval(b,xx,'identity',stats);

figure(1); clf; 
plot(xx,yfit,'r','linewidth',3); hold on; plot(x,yy,'k*'); axis square; 

%% For each 2 degrees increment find 25%, 75% percentile
for i=1:length(xx); 
   x0=max(xx(i)-1,min(xx)); x1=min(xx(i)+1,max(xx)); 
   ind=find(x>=x0 & x<=x1); 
   yu(i)=quantile(yy(ind),.95); 
   yd(i)=quantile(yy(ind),.05); 
end

axis tight
title('Temperature vs. BSP','fontname','arial','fontsize',13);; 
xlabel('Temperature [\circC]','fontname','arial','fontsize',12); 
ylabel('log[BSP/(1-BSP)]','fontname','arial','fontsize',12); 
set(gca,'fontname','arial','fontsize',11); 
%str=sprintf('y=8.57-0.335x\nCorrelation coeff: 0.98');
%text(30,3.7,str,'fontname','arial','fontsize',10); 
set(gcf,'color','w'); 

figure(2); clf; 
plot(xx,yu-yd,'m','linewidth',2); 
% fit gaussian to this
f=yu-yd; 
sig=linspace(0,20,1000); 
bestSig=inf; 
bestparms=[ 27.4588    4.5701    7.3348];
best=inf; 
for i=1:100000
    m=bestparms(1)+randn*.1; sig=abs(bestparms(2)+randn*.1); h=abs(bestparms(3)+randn*.1); 
    g=h*exp(-((xx-m)/2/sig).^2); 
    E=sum((f-g).^2); 
    if E<best
        best=E; 
        bestparms=[m sig h]; 
        gBest=g; 
        figure(2); clf; plot(xx,f,'k.',xx,g); drawnow
    end
end

%% Give rsquared value for fit
fb=mean(f); 
sst=sum((f-fb).^2);
ssr=sum((f-gBest).^2);
Rsq=1-(ssr/sst); 
disp(Rsq); 

figure(2); clf; 
plot(xx,f,'k*'); 
hold on; plot(xx,gBest,'r--','linewidth',3); drawnow
axis square
axis tight
title('Gaussian Fit to 95% Interquantile Ranges','fontname','arial','fontsize',13);; 
xlabel('Temperature [\circC]','fontname','arial','fontsize',12); 
ylabel('95% IQR of log[BSP/(1-BSP)]','fontname','arial','fontsize',12); 
set(gca,'fontname','arial','fontsize',11); 
set(gcf,'color','w'); 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng -r300 FigBSPIQR.png



%% Put fit onto figure 1
m=bestparms(1); sig=abs(bestparms(2)); h=abs(bestparms(3)); 
g=h*exp(-((xx-m)/2/sig).^2); 
figure(1); 
yu=yfit'+g/2; 
yd=yfit'-g/2; 
hold on; plot(xx,yu,'r--',xx,yd,'r--','linewidth',2); 
set(gcf, 'PaperPositionMode', 'auto');
print -dpng -r300 FigBSPvsTemp.png
