function [uevents,ueventtimes,ei,ARTINDS,BURSTINDS,BURST]=fcn_SegmentEEG(Fs,t,en,th,uevents,ueventtimes); 

%% Onset / offset: require change > th or < th that persists >0.5 sec
UTh=(Fs*0.5); DTh=Fs*.5; thE=th/2;
InBurst=0; j0=1; ctt=0; U=0; D=0; Burst=[];
for j=1:length(en); 
    % Actions when above threshold
    if (en(j)>th & ~InBurst & U==0); U=U+1; j0=j; end % just went above threshold -- possible start
    if (en(j)>th & ~InBurst & U>0); U=U+1; end % still above thresh
    if (en(j)>th & U>UTh); InBurst=1; end % have been above thresh >0.5 sec

    % Actions when below threshold
    if (en(j)<=thE & InBurst & D==0); D=D+1; j1=j; end % just dropped below -- possible end
    if (en(j)<=thE & InBurst & D>0); D=D+1; end % still below thresh
    if (en(j)<=thE & InBurst & D>DTh); ctt=ctt+1; Burst(ctt,:)=[j0 j1]; InBurst=0; end % stayed below for enough time to mark end
    
    if (en(j)<=th); U=0; end
    if (en(j)>th); D=0; end 
end

ind=find(uevents~=6); uevents(ind)=[]; ueventtimes(ind,:)=[]; %% clear old uevents and ueventtimes except for artifacts
BURST{1}=fcn_BurstOnsetOffsetTimes2(Fs,en,th); bb=BURST{1}; bb(end,2)=length(en); bb(1,1)=1; 

%% Put bursts into vector, then artifacts, all on top of suppressions. then read out into event log
ti=1:length(en); 
% 1 = burst, 2 = suppression, 6 = artifact
ei=2*ones(size(ti)); 

% put in bursts
for j=1:size(bb,1); 
   ind=bb(j,1):bb(j,2); 
   ei(ind)=1; 
end

% put artifacts on top
for j=1:length(uevents); 
   ind=find(t>=ueventtimes(j,1) &  t<=ueventtimes(j,2)); 
   ei(ind)=6; 
end

% create event log
uevents=[]; ueventtimes=[];
ct=0; j0=1; 
for j=2:length(en)
   if ei(j)~=ei(j-1); 
       j1=j-1; 
       ct=ct+1; 
       uevents(ct)=ei(j-1); 
       ueventtimes(ct,:)=[t(j0) t(j1)];
       j0=j; 
   end
end

%% Put event list also in index form-- referenced to eeg
ct=0; BURSTINDS=[];
for i=1:length(uevents)
   if uevents(i)==1; 
      ind=find(t>=ueventtimes(i,1) & t<=ueventtimes(i,2)); 
      ct=ct+1; 
      BURSTINDS(ct,:)=[ind(1) ind(end)];
   end
end

ct=0; ARTINDS=[];
for i=1:length(uevents)
   if uevents(i)==6; 
      ind=find(t>=ueventtimes(i,1) & t<=ueventtimes(i,2)); 
      ct=ct+1; 
      ARTINDS(ct,:)=[ind(1) ind(end)];
   end
end