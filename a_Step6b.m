clear all; clc; format compact;

%% Plot individual spectra for figures
td = '../Hypothermia-EEG-Data/'; 
pd = pwd; 

fileNo=7; 
FileNo=[1 2 3 4 5 7 8 10 12 13 14];
Tstart(FileNo)=[0   105 0   105 5   35  60  20  60  0   60]*60; % in sec
Tend(FileNo)=  [500 300 260 500 150 250 500 180 190 200 140]*60; % in sec
Tstart=Tstart(fileNo); Tend=Tend(fileNo); 
    

ct=0; SPECTRA=[]; S=[]; ct=0; xind=[];M=[];
cc=0; ch=0; 

for fileNo=FileNo(1:end)
    if ~(fileNo==7)
    %% Get detrended EEG and temperature
    % load(sprintf('INDIVBURSTDATA%d',fileNo)); % t s Fs en th smallTh BURST
    cd(td); 
    load(sprintf('CURATEDDATA%d',fileNo)); % t Fs s z tz ei en th t_temp Temp tbsp bsp ARTINDS BURSTINDS
    cd(pd); 

    Temp=(Temp-32)*5/9;
    %% Look at each identified burst, accept or reject [do only in first channel
    Nburst=size(BURSTINDS,1); 
    for j=2:Nburst; 
        t0=t(BURSTINDS(j,1)); t1=t(BURSTINDS(j,2)); 
        tm=(t0+t1)/2;
        
        dt=t(2)-t(1);
        Len=(BURSTINDS(j,2)-BURSTINDS(j,1))*dt;
        % Get temp
        [~,jj]=min(abs(t_temp-tm)); 
        Tt=Temp(jj);
%         disp([j Nburst Tt]); 

        i0=find(t>=max(0,(tm-2.5))); i0=i0(1); i1=find(t<min(max(t),tm+2.5)); i1=i1(end);     
        ss=s(i0:i1); tt=t(i0:i1)/60; 
        i0=max(1,i0-1500); i1=min(length(s),i1+1500); 
        sb=s(i0:i1); tb=t(i0:i1)/60; 
        
        movingwin=[2.68 .1]; 
        T=movingwin(1); K=3; 
        W=2*(K+1)/2/T;
        
        
        TW=T*W; K=2*TW-1; params.tapers=[TW K]; params.fpass=[.5 30]; params.pad=0; params.Fs=Fs; 
        [spect,tS0,f]=mtspecgram_mbw(ss,params,movingwin,0); % spectrum: Nt x Nf
        ct=ct+1; xind(ct)=size(S,2);
%         disp([T W K]);


            if (Len<1)
                disp(Len)
                cc=cc+1; 
                m=median(spect); %m=m/sum(m); 
                m2=mean(spect); %m2=m2/sum(m2); 
                figure(1); clf; 

                df=f(2)-f(1); ff=(-100:100)*df; 
                g=exp(-1/2*(ff/.251).^2); g=g/sum(g); 
                sss=conv(m,g,'same'); 
                subplot(131); imagesc(pow2db(spect')); axis xy
                subplot(132); plot(f,pow2db(m),f,pow2db(m2),'m',f,pow2db(sss),'k'); 
                subplot(133); plot(tb,sb); 
                S1(:,cc)=sss; 
%                 g=input('ok')                
            end
            if (Len>1 & Len<5)
                disp(Len)
                ch=ch+1; 
                m=median(spect); %m=m/sum(m); 
                m2=mean(spect); %m2=m2/sum(m2); 
                figure(1); clf; 

                df=f(2)-f(1); ff=(-100:100)*df; 
                g=exp(-1/2*(ff/.251).^2); g=g/sum(g); 
                sss=conv(m,g,'same'); 

                subplot(131); imagesc(pow2db(spect')); axis xy
                subplot(132); plot(f,pow2db(m),f,pow2db(m2),'m',f,pow2db(sss),'k'); 
                subplot(133); plot(tb,sb); 
                S2(:,ch)=sss;
%                 g=input('ok') 
            end
            
        end

    end    
    %%
end

m1=mean(S1'); 
m2=mean(S2'); 

cd(td); 
save(sprintf('ShortLong%d',100),'m1','m2','S1','S2'); 
cd(pd); 
