clear all; clc; format compact; 

%% Find differences between spectra for cold vs warm
clear all; clc; format compact; 
load FFF
% depends on: Step5A, 5B, 6 (run first)

%% Find differences between warm and cold BS spectra
% load WarmColdSpectra2 % m1 m2 S1 S2
fileNo=1; 

for fileNo=[1 3 7 100] %[1 2 3 4 5 7 8 10 12 13 14];
 
    load(sprintf('WarmCold%dS',fileNo)); 
    % remove outliers
    [S1,m1]=fcnRemoveOutliers(S1,m1);
    [S2,m2]=fcnRemoveOutliers(S2,m2);
    
    p=0.05; 
    % Get multiple median spectra for each patient GA
    Serr=fcn_ErrorBarsMultipleSpectra(m1,S1,p); 
    x=f; y1= pow2db(m1); 
    errBar1=[(y1- pow2db(Serr(1,:))); ( pow2db(Serr(2,:))-y1)]; 

    % Get multiple median power spectra for each patient -- Emergence
    Serr=fcn_ErrorBarsMultipleSpectra(m2,S2,p); %% stopped here
    x=f; y2= pow2db(m2); 
    errBar2=[(y2- pow2db(Serr(1,:))); ( pow2db(Serr(2,:))-y2)]; 

    %% Look at all traces -- visually assess variability
    figure(11); subplot(211); 
    % plot(f,pow2db(S1),'k'); subplot(212); plot(f,pow2db(S2),'k')
    plot(f,S1,'k'); subplot(212); plot(f,S2,'k')

    %% Do the two group test
    [dz,vdz,Adz,pval]=fcn_TGT_Spect(S1,S2,.01,'y',f);

    % Find rejected parts
    df=abs(f(2)-f(1)); 
    T=1.5; K=3; 
    W=2*(K+1)/2/T; 
    Nr= ceil(2*W/df); % number of consecutive pointwise rejections required to reject
    reject1=fcn_FindRejectedSegmentsB(Adz,Nr); 
    reject2=fliplr(fcn_FindRejectedSegmentsB(fliplr(Adz),Nr)); 
    reject=min((reject1+reject2)',1)'; 
    reject=~reject;
    fr=f; 

    %% Plots
    figure(1); clf; 
    plot(f,y1,'b',f,y2,'r'); hold on;
    z=nan(size(reject)); ind=find(reject==0); z(ind)=10; 
    plot(f,z,'k+'); 

    figure(10); clf; col='br'; 
    lineProps='b'; 
    plot(x,y1,'b','linewidth',2); hold on;
    H=shadedErrorBar(x,y1,(errBar1),lineProps,1);
    yu=y1+errBar1(1,:); yd=y1-errBar1(2,:); 
    % plot(f,yu,'b--',f,yd,'b--'); 

    lineProps='r'; 
    plot(x,y2,col(2),'linewidth',2);   
    H=shadedErrorBar(x,y2,(errBar2),lineProps,1);
    yu=y2+errBar2(1,:); yd=y2-errBar2(2,:); 
    % plot(f,yu,'r--',f,yd,'r--'); 

    set(gcf,'color',[1 1 1]);

    set(gca,'fontsize',14,'fontname','arial'); 
    set(gca, 'LineWidth', 1.2);
    xlabel('Frequency (Hz)','fontsize',14,'fontname','arial'); 
    ylabel('Power (dB)','fontsize',14,'fontname','arial');
    set(gca,'tickdir','out');
    box off

    accept=~reject; 
    yy=nan(size(accept)); 
    ind=find(accept==1); yy(ind)=nan; 
    ind=find(accept==0); yy(ind)=1; 
    fy=fr; ind=find(fr>0); fy=fy(ind); yy=yy(ind); 
    plot(fy,-60*yy,'k-','linewidth',5); 

    h=text(21,17,'15-22\circC'); set(h,'fontname','arial','fontsize',13);
    set(h,'color','b')
    h=text(21,10,'27-34\circC'); set(h,'fontname','arial','fontsize',13)
    set(h,'color','r')
    
    axis([.5 30 -65 20]); 
    set(gca,'xtick',[0.5 5 10 15 20 25 30]); 

    set(gcf, 'PaperPositionMode', 'auto');
    eval(['print -dpng -r300 Fig_WarmColdS' num2str(fileNo) '.png']);
    
    %**********************

    % axis tight
    figure(1); close
    figure(2); close
    figure(11); close   
%      g=input('ok'); 
end