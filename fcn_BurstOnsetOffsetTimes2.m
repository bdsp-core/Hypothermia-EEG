function Burst=fcn_BurstOnsetOffsetTimes2(Fs,en,th); 

%% Onset / offset: require change > th or < th that persists >0.5 sec
UTh=(Fs*0.25); DTh=Fs*.25; thE=th/2;
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