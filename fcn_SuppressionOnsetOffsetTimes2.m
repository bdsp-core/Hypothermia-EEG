function Supp=fcn_SuppressionOnsetOffsetTimes2(Fs,en,th); 

%% Onset / offset: require change > th or < th that persists >0.5 sec
dTh=(Fs*0.35); uTh=Fs*.35; 
InSupp=0; j0=1; ctt=0; d=0; u=0; Supp=[];
for j=1:length(en); 
    % Actions when below threshold
    if (en(j)<th(j) & ~InSupp & d==0); d=d+1; j0=j; end % just went below threshold -- possible start
    if (en(j)<th(j) & ~InSupp & d>0); d=d+1; end % still below thresh
    if (en(j)<th(j) & d>dTh); InSupp=1; end % have been below thresh >0.5 sec

    % Actions when above threshold
    if (en(j)>=th(j) & InSupp & u==0); u=u+1; j1=j; end % just rose above -- possible end
    if (en(j)>=th(j) & InSupp & u>0); u=u+1; end % still above thresh
    if (en(j)>=th(j) & InSupp & u>uTh); ctt=ctt+1; Supp(ctt,:)=[j0 j1]; InSupp=0; end % stayed above for enough time to mark end
    
    if (en(j)>=th(j)); d=0; end
    if (en(j)<th(j)); u=0; end 
end