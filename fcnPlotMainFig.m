function fcnPlotMainFig(fileNo)

%% Make main figures for paper -- do this only after step 11 is run
figure(2); clf

% tdS='..\..\BurstSuppressionPaper_Material\ForMP\Hypothermia';
% tdD='..\..\BurstSuppressionPaper_Material\ForMP\Hypothermia';
% tsc='..\..\BurstSuppressionPaper_Material\ForMP\Score BSup';
pd=pwd; 
td = '../Hypothermia-EEG-Data/'

FileNo=[1 2 3 4 5 7 8 10 12 13 14];
Tstart(FileNo)=[0   105 0   105 5   35  60  20  60  0   60]; % in sec
Tend(FileNo)=  [500 300 260 500 150 250 500 180 190 200 140]; % in sec


Tstart=Tstart(fileNo); 
Tend=Tend(fileNo); 

% for fileNo=[1 2 3 4 5 7 8 10 12 13 14]
    
    %% Load all data
    cd(td); 
    load(sprintf('spect%d',fileNo)); 
    
    S=spect{1}; %% Spectrogram

    eval(['load CURATEDDATA' num2str(fileNo)]); % t Fs s z tz ei th t_temp Temp tbsp bsp ARTINDS BURSTINDS
    
    cd(pd); 


    %% Zero times consistently for eeg, spectrogram; note: already done for Temp, bsp, infusions
    tspect=stimes/60; % spectroram time
    teeg=t/60; % eeg time
    tTemp=t_temp/60; 
    tBsp=tbsp/60; 
    Tend=min(Tend,max(teeg)); 
    tz=tz/60; 

%     ind=find(teeg>=Tstart(fileNo) & teeg<=Tend(fileNo)); s=s(ind); teeg=teeg(ind); 
%     ind=find(tspect>=Tstart(fileNo) & tspect<=Tend(fileNo)); S=S(:,ind); tspect=tspect(ind); 
    
    %% **********************************
    L=.1; W=.8; 
    
    %% **********************************
    % EEG        
    B=.85; H=.12;
    figure(2); 
    h=subplot('position',[L B W H]); set(h,'xtick',[]);  
    h=plot(teeg,s,'k'); 
    box off; 
    set(gca,'xtick',[]); 
    f=figure(2); 
    set(gca,'xcolor','w','tickdir','out'); 
    ylabel('Voltage [ \muV ]','FontName','Arial','FontSize',10);
    axis([Tstart Tend -100 100]);  
    set(gca,'FontName','Arial','FontSize',9); 

    %% **********************************
    % Binary data
    B=.82; H=.015;
    h=subplot('position',[L B W H]); set(h,'xtick',[]); set(h,'ytick',[]);  
    %ee=nan(size(ei)); 
    EI=ei; 
    ee=nan(size(EI));
    ind=find(ei==1); ee(ind)=-1; ind=find(ei==2); ee(ind)=1; ind=find(ei==6); ee(ind)=-.2;
    imagesc(tz,ones(size(tz)),ee,[-1 1]); colormap hot
    freezeColors
    axis([Tstart Tend 0 1]);  
    axis off 
    set(gca,'FontName','Arial','FontSize',9); 

    %% **********************************

    %% **********************************
    % Temp 
    xi=tTemp; yi=(Temp-32)/1.8;  %(T-32)/1.8;
    
    %% Add in infusion data
    load EEGSTARTTIME  % EEGSTARTTIME_ABS

    B=.1; H=.25;
    h=subplot('position',[L B W H]); 
    plot(xi,yi,'k','linewidth',2);
    ylabel('Temperature [ \circC ]','FontName','Arial','FontSize',12,'color','k'); 
    xlabel('Elapsed time [minutes]','FontName','Arial','FontSize',12,'color','k'); 
    axis([Tstart Tend 16 36]);  
    
    set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'ycolor'      , 'k'       ,...
    'LineWidth'   , 1         );
    set(gcf,'color',[1 1 1]); 
    grid on; %grid minor
    set(gca,'ytick',[16:5:34])
    set(gca,'FontName','Arial','FontSize',10); 

    %% **********************************

    % Spectrogram
    B=.65; H=.15;
    figure(2); 
    f=sfreqs; clim=[-20 10];
    h=subplot('position',[L B W H]); set(h,'xtick',[]);  
    colormap jet
    imagesc(tspect,f,pow2db(S),clim); axis xy 
    box off; 
    set(gca,'xtick',[]); 
    set(gca,'FontName','Arial','FontSize',10); 
    f=figure(2); 
    set(gca,'xcolor','w','tickdir','out'); 
    ylabel('Power [ dB ]','FontName','Arial','FontSize',10);
    axis([Tstart Tend min(sfreqs) 20]);
    set(gca,'FontName','Arial','FontSize',10); 

    %% **********************************

    %% **********************************
    % BSP
    B=.38; H=.25;
    h=subplot('position',[L B W H]); set(h,'xtick',[]); 
    plot(tBsp,bsp,'k','linewidth',2); 
    ylabel('BSP','FontName','Arial','FontSize',10); 
    axis([Tstart Tend -.01 1.2]);  
    set(gca,'xtick',[]); 
    set(gca,'xcolor','w','tickdir','out'); 
    grid on; 
    set(gca,'xtick',[0:0.5:10]); 
    set(gca,'ytick',[0 .5 1])
    set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YTick'       , 0:.25:1,...
    'LineWidth'   , 1         );
    set(gcf,'color',[1 1 1]); 
    set(gca,'FontName','Arial','FontSize',10); 
    
    %% ***********************************
    % Mark times of ECS
    tt=t/60; dt=tt(2)-tt(1);     
    ind=find(tt>Tstart & tt<Tend); tt=tt(ind); 
    Ei=ei(1:length(tt)); ct=0; S=s(ind); En=en(ind);
    A=zeros(1,length(S)); ECI=nan(size(A));
    for j=1:length(S); 
        if (En(j)<2 & ~(Ei(j)==6)); 
           ct=ct+1; A(j)=dt*ct;
        elseif (En(j)>=2 & ~(Ei(j)==6)); 
            ct=0; 
        else 
            % do nothing
        end
        if A(j)>3; ECI(j)=1; end
    end
    hold on; plot(tt,ECI,'r-','linewidth',3)
    %% **********************************
    
    %% Save figure

    FileName=['FigEEG_BSR_Temp_' num2str(fileNo)]; 
    set(gcf, 'PaperPositionMode', 'auto'); 
    print('-dpng','-r300',FileName);
    %% **********************************

