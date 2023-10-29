function [SPECTRA,BurstWaveform,tb,bl,bt,f]=fcnGetIsolatedBursts(t,s,Fs,en,th,t_temp,TEMP,smallTh,BURST)

%% Look at each identified burst, accept or reject [do only in first channel
Burst=BURST{1};
Nburst=size(Burst,1); 
ct=0; SPECTRA=[]; ta=[]; Temp=[]; bt=[]; bl=[]; pp=[]; sss=[];
for j=1:Nburst; 
    t0=t(Burst(j,1)); t1=t(Burst(j,2)); 
    i0=find(t>=max(0,(t0-20))); i0=i0(1); i1=find(t<min(max(t),t1+20)); i1=i1(end);    
    t0i=t(i0); t1i=t(i1); 
    ss=s(i0:i1); tt=t(i0:i1); ens=en(i0:i1); 
    xx=[t0 t1]; yy=[th th]; 
    figure(1); clf; subplot(2,3,1:2); plot(tt,-ss,'k',tt,ens,'r',xx,yy,'b*',tt,ones(size(tt))*th,'m--'); 
    axis([min(tt) max(tt) -25 25]); 
    set(gcf,'color','white'); box off; xlabel('seconds'); ylabel('uv'); 

    % extend edges further
    aa=smallTh; 
    idx=1:Burst(j,1); % R side of burst
    j0=find(en(idx)<aa);
    idx=Burst(j,2):length(en); 
    j1=find(en(idx)<aa); % L side of burst
    if ~isempty(j0) & ~isempty(j1); 
        j0=j0(end); j1=j1(1)+Burst(j,2); Fs=round(Fs); 
        j0=max(1,j0-Fs); j1=min(length(t),j1+Fs); 
    else
        j0=Burst(j,1); j1=Burst(j,2); 
    end
    xx=[t(j0) t(j1)]; yy=[th th];
    subplot(2,3,1:2); hold on; plot(xx,yy,'go'); 

    % plot extracted burst
    subplot(2,3,4:5); plot(t(j0:j1),-s(j0:j1),'k');
    axis([min(tt) max(tt) -25 25]); 
    set(gcf,'color','white'); box off; xlabel('seconds'); ylabel('uv');     

    BurstWaveform{j}=-s(j0:j1); 

    movingwin=[4 .1]; 
    T=movingwin(1); 
    W=1/2; 
    TW=T*W; 
    K=2*TW-1; %disp([T W K]); 
    params.tapers=[TW K]; 
    params.fpass=[.5 30]; 
    params.pad=0; 
    params.Fs=Fs; 

    [S,tS0,f]=mtspecgram_mbw(ss,params,movingwin,0); % spectrum: Nt x Nf
    subplot(233); imagesc(tS0,f,pow2db(S')); axis xy;

    pp{j}=S; 
    
    % get spectrum of burst
    tf0=t0-t0i; tf1=t1-t0i; 
    ind=find(tS0>=tf0 & tS0<=tf1);     
    sf=mean(S(ind,:)); 

    try
        subplot(236); plot(f,pow2db(sf)); axis xy;
    catch
    end
    drawnow

    % collect spectra & times; save temp with it
    ta=[ta (t0+t1)/2]; % time of burst relative to beginning of recording
    [~,jj]=min(abs(t_temp-ta(end))); Temp=[Temp TEMP(jj)];
    try
    SPECTRA=[SPECTRA; sf];    
    catch
    keyboard
    end
    bl=[bl (t1-t0)];
    bt=[bt TEMP(jj)];
end    
tb=ta; 
