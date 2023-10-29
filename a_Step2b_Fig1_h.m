clear all; clc; format compact; 

%% Compute spectra for each temperature bin from multiple bursts of different lengths
% only use those >4 sec duration
cn=0;     
ff=figure(4); clf;
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.1 0.01], [0.1 0.01]);
td = '../Hypothermia-EEG-Data/'; 

for CaseNo=3 %[1 2 3 4 5 6 7 8 10 12 13 14 ]; %
    clear bt pp bl f BURST TBURST

    filename = ['BurstSpectraTempIndiv_' num2str(CaseNo)];
    full_path = [td filename];
    load_cmd = ['load ''' full_path ''''];
    eval(load_cmd); % ' pp bt bl f BURST TBURST']
   
    bt=(bt-32)*5/9; % convert to celsius
    for j=1:length(bl); 
        bl(j)=(1/200)*length(BURST{j}); 
        tb=mean(TBURST{j}); 
    end
    Tburst=bt; 
    Burst=BURST;
    %% Reduce number
    ct=0; 
    %Burst=Burst(60:end); Tburst=Tburst(60:end); 
    for i=1:length(Burst); 
        if ~mod(i,2)
            ct=ct+1; 
            B{ct}=Burst{i}; 
            Tb(ct)=Tburst(i); 
        end
    end
    Burst=B; Tburst=Tb; 
    %% Plot bursts on one page
    Fs=200; nsec=45; Nts=Fs*nsec; % max seconds per row
    i0=0; y=50; figure(1); clf;

    % Index colors to temps
    col=colormap(jet); col=flipud(col); tmx=32; tmn=18; 
    Tnums=linspace(tmx,tmn,size(col,1)); 
    for i=1:length(Tburst); [~,TburstCol(i)]=min(abs(Tnums-Tburst(i))); end

    set(gcf,'color',[0 0 0]); 
    for i=1:length(Burst); 
       f=Burst{i}; %if length(f)>1000; f=f(1:1000); end
       if length(f)>0; 
       f=f-mean(f); %f=f/std(f); 
       i1=i0+length(f)-1; 
       if i0>Nts; i0=0; i1=i0+length(f)-1;  y=y-40; end   
       x=i0:i1; f=f+y; 
       figure(1); hold on; plot(x,f,'k','color',col(TburstCol(i),:),'linewidth',1); axis off; 
       %text(x(1),f(1)+5,sprintf('%.1f',Tburst(i)),'color','white','fontsize',8); 
       i0=i1+100;  
       end
    end
    
%     g=input('ok')
    
    
end

%%
set(ff,'color',.5*[1 1 1])
figure(1); set(gcf,'color',0*[1 1 1]); 
h=colorbar; 
ax=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3)*0.95 pos(4)]);
%%
pos=get(gca,'pos');
hc=colorbar('location','east','position',[.75 .15 .02*pos(3) .85*pos(4)]);
set(hc,'xaxisloc','bottom');
set(hc,'yaxisloc','right');
set(hc,'box','off')
% 
% set(gcf, 'PaperPositionMode', 'auto');
% print -dpng -r300 Fig1_Spectra.png
caxis([tmn tmx])
set(get(hc,'xlabel'),'string','Temp (\circ C)','fontname','arial','fontsize',12,'rotation',0);
xx=(Nts-10*Fs):Nts;

%%
yy=(-1570)*ones(size(xx)); 
plot(xx,yy,'w');
xx=[min(xx) min(xx)];
yy=[min(yy) (min(yy)+50)];
plot(xx,yy,'w');
text(6400,-1550,'\mu v','color','w','fontname','arial','fontsize',12);
text(7000,-1595,'10 sec','color','w','fontname','arial','fontsize',12);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'InvertHardCopy', 'off');
print -dpng -r300 Fig3_ColorPlots.png
