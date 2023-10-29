clear all; clc; format compact; 

%% Compute burst length distributions for cardiac surgery hypothermia data
% Plots: (1) overall distribution; (2) distribution for hot vs cold temps

td = '../Hypothermia-EEG-Data/'; 

ct=0; dt=1/200; 
for CaseNo=[1 2 3 4 5 6 7 8 10 12 13 14 ]; % [4 5 6 7 8 10 12 13]
    disp(['Working on record ' num2str(CaseNo)]); 
    full_path = [td 'SPECTRA' num2str(CaseNo)];
    load(full_path); 
    Nb=length(bt);    
    for j=1:Nb
        L=bl(j); TT=bt(j); 
        d=BurstWaveform{j}; 
        a=max(d)-min(d); 
        if ~isempty(d); 
%            d=fcn_CleanUp_d(d); 
            if (L>0 & L<30 & ~isnan(TT) & a<300); 
                ct=ct+1; Len(ct)=L; T(ct)=TT; A(ct)=max(d)-min(d); 
            end
        end
    end
end
T=(T-32)/1.8; 

%% Make cdfs for some temperature devisions
figure(1); clf; 
bins=[17 22; 27 32];
for i=1:size(bins,1);
   ind=find(T>=bins(i,1) & T<bins(i,2)); 
   cL{i}=Len(ind); 
   cA{i}=A(ind);
end

%% Plot cdfs
figure(2); clf; 
col='br';
for i=1:size(bins,1); 
   figure(1); 
   [h,stats,xA{i},yA{i}]=cdfplot_mbw(cA{i}); 
   
   figure(2); hold on;  
   x=xA{i}; y=yA{i}; 
   x(end)=250; y(end)=1; 
   plot(x,y,col(i),'linewidth',2); hold on
   axis([0 250 0 1.05]); 
   xlabel('Burst amplitude (microvolts)','fontsize',14,'fontname','arial');
   ylabel('CDF','fontsize',14,'fontname','arial');
   set(gca,'fontsize',14,'fontname','arial');
   set(gca,'tickdir','out','XMinorTick','on','YMinorTick','on'); 
   box off
end
figure(2); set(gcf,'color',[1 1 1]); box off; 
h=legend('17-22 \circC','27-32 \circC'); set(h,'location','southeast'); set(h,'box','off')

%% Plot D from KS test

% For amplitudes
x1=cA{1}; x2=cA{2}; 
[H,P,KSSTAT]=kstest2(x1,x2,0.05,'unequal'); 
x=sort([xA{1}; xA{2}]); mx=-inf; 
N1=length(cA{1}); 
N2=length(cA{2}); 
for i=1:length(x); 
    c1=sum(cA{1}<=x(i))/N1; c2=sum(cA{2}<=x(i))/N2; 
    d=abs(c1-c2); 
    if d>mx; xmx=x(i); mx=d; y0=min(c1,c2); y1=max(c1,c2); end
end
figure(2);
xx=[xmx xmx]; yy=[y0 y1]; plot(xx,yy,'k--',xx,yy,'ko','markersize',5','markerfacecolor','k'); 
[xx,yy] = ds2nfu(xx(1)+.5,mean(yy-.1));
annotation('textbox',[xx(1) yy(1) .1 .1],'string',sprintf('KS test result: \n D = %.2f, p<10^{-97}',mx),...
    'fontname','arial','fontsize',14,...
    'backgroundcolor',[1 1 1],...
    'edgecolor',[1 1 1]);

%text(xmx+1,(y0+y1)/2,sprintf('D = %.2f',d),'fontname','arial','fontsize',14,'background',[1 1 1]); 

set(gcf, 'PaperPositionMode', 'auto');
print -dpng -r300 Fig2_LengthAmplitudeCDFs.png
