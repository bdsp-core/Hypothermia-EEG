function [tbsp,bsp,z]=fcnGetBSP(tz,ei,fileNo)

%% Compute binary sequence and bsp
%  z=fcnCreateBinaryVectorV2(ei);

z=nan(size(ei)); 
Nt=length(z); 

ind=find(ei==1); z(ind)=1; % put bursts in first
ind=find(ei==2); z(ind)=0; % put suppressions in next

% local imputation -- replace nans by majority in 1 minute neighborhood
for i=1:Nt
   if isnan(z(i)); 
      i0=max(1,i-6000); i1=min(Nt,i+6000); 
      if sum(z(i0:i1)==1)>sum(z(i0:i1)==0); z(i)=1; else z(i)=0; end
   end
end

if fileNo==5; 
        ind=find((tz/60)>=140); z(ind)=1; 
end
% Delete data before start time -- is garbage, and biases calculations
% ind=find(t<=Tstart(fileNo)); t(ind)=[]; z(ind)=[]; 
% Delete data after end time -- is garbage, and biases calculations

% ind=find(t>=Tend); t(ind)=[]; z(ind)=[]; 
zz=z(1:50:end); tz=tz(1:50:end); 
%zz=z; tz=t; 
dt=tz(2)-tz(1); 
Fs=1/dt; 

%% Run EM code
binomial_trials = 10; tstart = tic; %Fs=200; 
[Responses, ptile,xnew,signewsq,A,p, x, ss] = BSP_EM_estim(~zz,Fs,binomial_trials);
tz=linspace(min(tz),max(tz),length(ptile.pmode)); 
lower=ptile.p050; upper=ptile.p950;
bsp=ptile.pmode; 
tbsp=tz; 