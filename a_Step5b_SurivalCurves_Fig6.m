clear all; clc; format compact; 

% Make "survival curves"
td = '../Hypothermia-EEG-Data/'; 
pd = pwd; 

cd(td); 
load SurvivalData % Ev1 Ev0
cd(pd); 

nb=50; nL=100; 
T=linspace(20,35,nb); 
L = linspace(eps,60,nL);
S0=zeros(nL,nb); S1=S0; 
b=(T-15)/(35-15);

for i=1:nb; 
    %[h0,h1,wc0,wc1,wp0,wp1,a0(i),b0(i),a1(i),b1(i)]=fcn_GetWeibullFits(L,b(i),Ev1,Ev0);
    T0=(Ev0(:,2)+Ev0(:,3))/2; bT0=(T0-min(T))/(max(T)-min(T)); L0=Ev0(:,1); 
    T1=(Ev1(:,2)+Ev1(:,3))/2; bT1=(T1-min(T))/(max(T)-min(T)); L1=Ev1(:,1);
    t0=fcn_Get_FailTimes(bT0,b(i),L0); 
    t1=fcn_Get_FailTimes(bT1,b(i),L1); 
    [f0,x0]=ecdf(t0); 
    [f1,x1]=ecdf(t1); 
    x0=x0(2:end); x1=x1(2:end); f0=f0(2:end); f1=f1(2:end); 
    f0=interp1(x0,f0,L); f1=interp1(x1,f1,L); 
    ind=find(isnan(f0)); f0(ind)=0; 
    ind=find(isnan(f1)); f1(ind)=0; 
    % get smaller with increasing L
    for j=2:nL
        f0(j)=max(f0(j),f0(j-1)); 
        f1(j)=max(f1(j),f1(j-1)); 
    end
    
    S0(:,i)=f0; 
    S1(:,i)=f1; 
end

%% Figure
figure(1); clf; hold on; 
subplot(211); imagesc(b,L,(S0)); axis xy; axis([min(b) max(b) 0 20]); 
subplot(212); imagesc(b,L,(S1)); axis xy; axis([min(b) max(b) 0 20]); 

%% Prelim fit
%% Use fitted parameters to make plot

FitAgain=0; 

if FitAgain==1
%% Surface fits
% Bursts
parms0=[ 7.4462    1.3066    0.8201];
S=S1; [xx,fval] = fminsearch(@(parms) fcn_FitWeibull0(parms,S,L,b),parms0);
c1=abs(xx); 

% Suppressions
parms0=abs([  9.5970    1.2955    0.8686]); % 
S=S0; [xx,fval] = fminsearch(@(parms) fcn_FitWeibull1(parms,S,L,b),parms0);
c0=abs(xx); 


save WEIBULLSURFACEPARMS c0 c1
else
    load WEIBULLSURFACEPARMS
end


%% Figure for poster

c=c0; A=c(1)*exp(c(2)*b.^2); B=c(3); A=max(A,eps); B=max(B,eps); 
for i=1:length(b);    wc = wblcdf(L,A(i),B); s(:,i)=wc; end 
s0=s; 
m0=A*log(2)^(1/B); 

c=c1; A=c(1)*exp(c(2)*(1-b).^2); B=c(3); A=max(A,eps); B=max(B,eps); 
for i=1:length(b);    wc = wblcdf(L,A(i),B); s(:,i)=wc; end 
s1=s; 
m1=A*log(2)^(1/B); 



figure(1); clf; hold on; 
h=subaxis(2,2,1,'SpacingVert',0,'SpacingHoriz',0,'PR',0,'PL',.01'); 
    contourf(T,L,(S0),[0:.1:1.01]); axis xy; axis([min(b) max(b) 0 60]); 
    hold on; plot(T,m1,'k--','linewidth',3);
    set(h,'xtick',[]); 
    set(h,'fontname','arial','fontsize',15); 
    axis([min(T) max(T) 0 59]);
    set(gca,'tickdir','out'); 
    ylabel('Duration (sec)','fontname','arial','fontsize',20); 
    title('Suppression Duration CDF','fontname','arial','fontsize',20);

h=subaxis(2,2,2,'SpacingVert',0,'SpacingHoriz',0,'PR',.01); 
    contourf(T,L,(S1),[0:.1:1.01]); axis xy; axis([min(b) max(b) 0 60]); set(gca,'xtick',[]); set(gca,'ytick',[]); ; 
    hold on; plot(T,m0,'k--','linewidth',3);
    set(gca,'fontname','arial','fontsize',15); 
    axis([min(T) max(T) 0 59]);
    set(gca,'tickdir','out'); 
    title('Burst Duration CDF','fontname','arial','fontsize',20);
    set(h,'xtick',[]); 

h=subaxis(2,2,3,'SpacingVert',0,'SpacingHoriz',0,'PR',0,'PL',.01'); 
    contourf(T,L,(s1),[0:.1:1.01]); axis xy; axis([min(b) max(b) 0 60]);
    hold on; plot(T,m1,'k--','linewidth',3);
    set(gca,'fontname','arial','fontsize',15); 
    axis([min(T) max(T) 0 59]);
    set(gca,'xtick',[min(T):2:max(T)]);     
    set(gca,'tickdir','out'); 
    ylabel('Duration (sec)','fontname','arial','fontsize',20); 
    xlabel('Temperature [\circC]','fontname','arial','fontsize',20);
    
    subaxis(2,2,4,'SpacingVert',0,'SpacingHoriz',0,'PR',0.01,'PT',.01); 
    contourf(T,L,(s0),[0:.1:1.01]); axis xy; axis([min(b) max(b) 0 60]); 
    hold on; plot(T,m0,'k--','linewidth',3);
    set(gca,'fontname','arial','fontsize',15); 
    axis([min(T) max(T) 0 59]);
    set(gca,'xtick',[min(T):2:max(T)]);     
    set(gca,'tickdir','out'); 
    xlabel('Temperature [\circC]','fontname','arial','fontsize',20);
    set(gca,'ytick',[]); 
    
colormap hot
set(gcf,'color',[1 1 1]); 

FileName='FigCDFs'; 
set(gcf, 'PaperPositionMode', 'auto'); 
print('-dpng','-r300',FileName);

%%
figure(2); 
s0(1,:)=0; 
s0(100,:)=1.01; 
    contourf(T,L,(s0),[0:.01:1]); axis xy; axis([min(b) max(b) 0 60]);  
    set(gca,'fontname','arial','fontsize',15); 
    axis([min(T) max(T) 0 59]);
    set(gca,'xtick',[min(T):2:max(T)]);     
    set(gca,'tickdir','out'); 
    xlabel('Temperature [\circC]','fontname','arial','fontsize',20);
    set(gca,'ytick',[]); 

    %
h=colorbar; colormap hot
set(gcf,'color',[1 1 1]); 
set(h,'ytick',[0 .25 .5 .75 1],'yticklabel',{'0.0','0.25','0.5','0.75','1.0'}); 
set(h,'fontname','arial','fontsize',15);
FileName='FigCDFsColorBar'; 
set(gcf, 'PaperPositionMode', 'auto'); 
print('-dpng','-r300',FileName);

%% Get R squared values
St0=sum(sum((S0-mean(S0(:))).^2));
St1=sum(sum((S1-mean(S1(:))).^2));
Sr0=sum(sum((s0-S0).^2));
Sr1=sum(sum((s1-S1).^2));
Rs0=1-Sr0/St0; 
Rs1=1-Sr1/St1;
disp([Rs0 Rs1]); 

%% Show median burst duration at 20, 30; same for suppressions
disp('***********'); 
ind=find(T>20); ind=ind(1); 
disp([m0(ind) m1(ind)]); 
ind=find(T>30); ind=ind(1); 
disp([m0(ind) m1(ind)]); 


