clear all; clc; format compact;  

td = '../Hypothermia-EEG-Data/'; 

fileNo=[1 2 3 4 5 7 8 10 12 13 14];
Tstart(fileNo)=[0   105 0   105 5   35  60  20  60  0   60]*60; % in sec
Tend(fileNo)=  [500 300 260 500 150 250 500 180 190 200 140]*60; % in sec
y=[]; x=[]; lenS=[]; ct=0; len=[]; ct=0; ctt=0; Ev1=[]; Ev0=[]; 

for fileNo=[1 2 3 4 5 7 8 10 12 13 14] % [1 2 3 4 5 7 8 10 12 13 14]
    disp(fileNo);

    full_path = [td 'SPECTRA' num2str(fileNo)]; load(full_path); 
    full_path = [td 'CURATEDDATA' num2str(fileNo)]; load(full_path); 

%     eval(['load SPECTRA' num2str(fileNo)]); % get: SPECTRA BurstWaveform tb bl bt f']);
%     eval(['load CURATEDDATA' num2str(fileNo)]); % t Fs s z tz ei en th t_temp Temp tbsp bsp ARTINDS BURSTINDS BURST']); 
    cd=0; dt=tz(2)-tz(1);
    ti=interp1(t_temp,Temp,tz); 
    Ti=(ti-32)*5/9; 
    %% Get suppression lengths
    tt=0; 
    for j=2:length(ei); 
        if ei(j)==2 & tt==0; ; 
            tt=tt+dt; i0=j; 
        elseif ei(j)==2 & tt>0
            tt=tt+dt; 
        elseif (ei(j)~=2 & ei(j-1)==2)
            ct=ct+1; i1=j; 
            Ev0(ct,:)=[tt min(Ti(i0:i1)) max(Ti(i0:i1))];
            tt=0; 
        end        
    end        
    
    tt=0; 
    for j=2:length(ei); 
        if ei(j)==1 & tt==0; ; 
            tt=tt+dt; i0=j; 
        elseif ei(j)==1 & tt>0
            tt=tt+dt; 
        elseif (ei(j)~=1 & ei(j-1)==1)
            ctt=ctt+1; i1=j; 
            Ev1(ctt,:)=[tt min(Ti(i0:i1)) max(Ti(i0:i1))];
            tt=0; 
        end        
    end        
end

pd = pwd; 
clear cd
cd(td); 
save SurvivalData Ev1 Ev0
cd(pd);