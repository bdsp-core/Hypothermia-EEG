function a_FIG5_IndividualBurstsInThrees(fileNo); 

% clear all; clc; format compact;

%% Plot individual spectra for figures -- plot 3 next to each other at specified times
% fileNo=3; 
FileNo=[1 2 3 4 5 7 8 10 12 13 14];

Tstart(FileNo)=[0   105 0   105 5   35  60  20  60  0   60]; % in sec
Tend(FileNo)=  [500 300 260 500 150 250 500 180 190 200 140]; % in sec
Tstart=Tstart(fileNo); Tend=Tend(fileNo); 
Times(1,:)=[24 107 198];
Times(3,:)=[37 124 203];
Times(7,:)=[117 183 210];

LL=[.13 .4 .67];
WW=.25*ones(1,3);


%% Get detrended EEG and temperature
td = '../Hypothermia-EEG-Data/';
fileName = [td 'INDIVBURSTDATA' num2str(fileNo)];
load(fileName); 
% load(sprintf('INDIVBURSTDATA%d',fileNo)); % t s Fs en th smallTh BURST

colormap jet
%% Look at each identified burst, accept or reject [do only in first channel
Burst=BURST{1};
Nburst=size(Burst,1); 
ct=0; SPECTRA=[]; 
figure(1); clf; 
for k=1:3; 
    for j=1:Nburst;     
        t0=t(Burst(j,1)); t1=t(Burst(j,2)); 
        tm=(t0+t1)/2;
        i0=find(t>=max(0,(tm-30))); i0=i0(1); i1=find(t<min(max(t),tm+30)); i1=i1(end);     
        ss=s(i0:i1); tt=t(i0:i1)/60;         
        if (round(tt(1))==Times(fileNo,k)); 
            disp(5*(bt(j)-32)/9); 
            movingwin=[2 .1]; 
            T=movingwin(1); W=1; 
            TW=T*W; K=2*TW-1; params.tapers=[TW K]; params.fpass=[.5 30]; params.pad=0; params.Fs=Fs; 
            %disp([T W K]); 
            [S,tS0,f]=mtspecgram_mbw(ss,params,movingwin,0); % spectrum: Nt x Nf
            fcnMakeTheFig(LL(k),WW(k),tt,ss,tS0,f,S,k); 
        end
   end    
end    

FileName=sprintf('FigBurstExs_Case%d',fileNo); 
set(gcf, 'PaperPositionMode', 'auto'); 
print('-dpng','-r300',FileName);


%% Make figure
function fcnMakeTheFig(L,W,tt,ss,tS0,f,S,k)
figure(1); hold on  
    B=.58; H=.35; h=subplot('position',[L B W H]);   
    h=plot(tt,-ss,'k'); 
    axis([min(tt) max(tt) -50 50]); 
    set(gca,'xtick',[],'tickdir','out','xcolor','w','fontname','arial','fontsize',10); 
    set(gcf,'color','white'); box off; 
    if k==1; 
        ylabel('Voltage [\muV]','fontname','arial','fontsize',10); 
    else
        set(gca,'ycolor','w'); 
    end

    B=.22; H=.35; a=subplot('position',[L B W H]);
    tS0=linspace(tt(1),tt(end),length(tS0)); 
    imagesc(tS0,f,pow2db(S'),[-20 20]); axis xy;
    xlabel('Elapsed time [minutes]','fontname','arial','fontsize',10);
    if k==1; 
        ylabel('Power [dB]','fontname','arial','fontsize',10); 
    else
        set(gca,'ycolor','w');
        set(gca,'yticklabel',''); 
    end
    set(gca,'tickdir','out'); 
    set(a,'box','off','color','none','fontname','arial','fontsize',10);
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    axes(a); 
%     linkaxes([a b])
    t0=round(tt(1))+.25; t1=t0+.5; 
    xl={num2str(t0),num2str(t1)}; 
    set(gca,'xtick',[min(tt+.25) max(tt-.25)],'xticklabel',xl); 
    drawnow