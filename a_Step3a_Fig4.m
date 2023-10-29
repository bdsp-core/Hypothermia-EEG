clear all; clc; format compact;  

%% Creates event files, binary vector files, spectra files
td = '../Hypothermia-EEG-Data/'; 


% These specify the useful part of the data
% Note time convention: refer all to eeg elapsed time;
fileNo=[1 2 3 4 5 7 8 10 12 13 14];
Tstart(fileNo)=[0   105 0   105 5   35  60  20  60  0   60]*60; % in sec
Tend(fileNo)=  [500 300 260 500 150 250 500 180 190 200 140]*60; % in sec
figure(1); clf; 

N=200; 
xt=linspace(18,32,N); 
xb=linspace(0,1,N); 
[X,Y]=meshgrid(xt,xb);
sigt=1; sigb=.05; 
x=[]; y=[];
figure(1); clf;
I=zeros(N,N); ct=0; 
for fileNo=[1 2 3 4 5 7 8 10 12 13 14] %14
    disp(fileNo);b=[]; 

    full_path = [td sprintf('CURATEDDATA%d',fileNo)];
    load(full_path); % t Fs s z tz ei en th t_temp Temp tbsp bsp ARTINDS BURSTINDS
    
    Temp=(Temp-32)*5/9;
    for j=1:length(Temp); 
       [~,jj]=min(abs(t_temp(j)-tbsp) );
       b(j)=bsp(jj);
       if Temp(j)<32;
           for i=1:length(xb); 
               for k=1:length(xt);                
                   I(i,k)=I(i,k)+exp(-1/2*(((xb(i)-b(j))/sigb)^2+((xt(k)-Temp(j))/sigt)^2)); 
               end; 
           end       
       end
    end
    ct=ct+1; xx{ct}=Temp; yy{ct}=b; 
    x=[x Temp]; y=[y b];
end
    
for j=1:N
   h=I(:,j); h=h/sum(h); I(:,j)=h;  
end
figure(1); imagesc(xt,xb,I); axis xy; drawnow;

pd = pwd; 
cd(td); 
save BSPTEMP xx yy x y
cd(pd); 