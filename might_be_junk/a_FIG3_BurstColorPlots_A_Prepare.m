clear all; clc; format compact; 

%% Compute spectra of suppressions in cardiac surgery hypothermia data
td = '/Users/mbw/cdac Dropbox/Datasets/zz_Hypothermia/Hypothermia-EEG-Data/'; 

ct=0; PP=[]; PP0=[]; PP1=[]; BT=[]; BL=[]; % group ensemble    
load HPF

for CaseNo=[1 2 3 4 5 6 7 8 10 12 13 14 ]; % [4 5 6 7 8 10 12 13]
    pp=[]; pp0=[]; pp1=[]; bt=[]; bl=[]; ctt=0; TBURST=[]; BURST=[]; % individual cases
    disp(['Working on record ' num2str(CaseNo)]); 
    
    filename = ['Data_' num2str(CaseNo)];
    full_path = [td filename];
%     load_cmd = ['load ''' full_path ''''];
%     eval(full_path); % ' pp bt bl f BURST TBURST']

    load(full_path); 
    Nb=size(Data.data,2);    
    
    
    for j=1:Nb
        d=Data.data{j}; 
        if ~isempty(d) & length(d)<(200*20); 
            d=locdetrend(d,200); 
            d(1:10)=0; d(end-10:end)=0;
            %d=filtfilt(HPF,1,d);
            t=Data.time{j}; 
            TempBurst(j)=Data.Temp{j}; 
            d(1:10)=0; d(end-10:end)=0;
       
            %****************
            % Filtering
            Fs=200; dt=1/200; tt=-.1:dt:.1; 
            % remove very low frequency 
%             sig=.4; g=exp(-1/2/sig^2*tt.^2); g=g/sum(g); fs=conv(d,g,'same'); fb=d-fs; 
%             d=fb; 
            % remove high frequency noise
            sig=.01; g=exp(-1/2/sig^2*tt.^2); g=g/sum(g); fh=conv(d,g,'same'); fb=fh; 
            d=fb; 
            %****************
        
            %% Compute spectrogram
            TW=3; K=2*TW-1;  % [default 5,9; 5,3 is same as Laura Lewis' paper]
            params.tapers=[TW K]; params.fpass=[1 15]; params.pad=1; movingwin=[2 .2]; params.Fs=200; 
            [S,tS0,f]=mtspecgramc(d,movingwin,params); % spectrum: Nt x Nf
            disp(f(2)-f(1)); % Frequency resolution
            %% Get burst coordinates
            indB=Data.burstIdx{j}; tburst=t(indB); NoSec=max(tburst)-min(tburst); burst=d(indB); 
            % Burst length
            L=length(indB)/200;
                
            %% Extract part of spectrogram within the burst
            idx=find(tS0>=min(tburst)&tS0<=max(tburst)); Sb=S(idx,:); tS=tS0(idx); 
        
            %% Average over time
            nh=round(length(tS0)/2); ii=round(nh/2):round(length(tS0-nh/2));
            P=sum(S(ii,:)); nt=round(size(Sb,1)/2); 

            % "Whiten" -- divide by 1/f spectrum
            
            pf=1./(f+.0001).^2; 
            pf=1; 
            P=P./pf; P=P/sum(P); 

            %% Plots
%             figure(1); set(gcf,'color',[1 1 1]); 
%             subplot(221); imagesc(tS0,f,log(S'));axis xy; axis([min(t) max(t) 3 15]); axis off; 
%             subplot(223); plot(t,d,'k',tburst,burst,'r'); axis([min(t) max(t) -50 50]); xlabel('t [sec]'); ylabel('uv'); box off;  
%             subplot(222); semilogy(f,P,'k'); %,f,P0,'r',f,P1,'b'); 
%             axis([min(f) max(f) 0 .05]);
%             box off; set(gcf,'color',[1 1 1]); 
%             xlabel('f [Hz]'); ylabel('log(P)'); 
%             drawnow; 
%             disp(TempBurst(j))
            %g=input('ok');  

            %% Store burst, temp for further analysis
            ct=ct+1; PP(ct,:)=P; BT(ct)=TempBurst(j); BL(ct)=L;  % group
            ctt=ctt+1; pp(ctt,:)=P;bt(ctt)=TempBurst(j); bl(ctt)=L; 
            BURST{ctt}=burst'; TBURST{ctt}=tburst;
            % individual cases
            %g=input('ok');  
        end
    end
    
    %% Save data for this patient
%     str=['save BurstSpectraTempIndiv_' num2str(CaseNo) ' pp bt bl f BURST TBURST']; 
%     eval(str); 
    full_path = [td 'BurstSpectraTempIndiv_' num2str(CaseNo)];
    save(full_path, 'pp', 'bt', 'bl', 'f', 'BURST', 'TBURST')

end

% Save aggregated data
full_path = [td 'BurstSpectraTempAll'];
save(full_path, 'PP', 'BT', 'BL', 'f')
% str=['save BurstSpectraTempAll' ' PP BT BL f']; eval(str); 
