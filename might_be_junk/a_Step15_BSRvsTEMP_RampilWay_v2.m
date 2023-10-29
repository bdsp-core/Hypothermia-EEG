clear all; clc; format compact; 
td = '/Users/mbw/cdac Dropbox/Datasets/zz_Hypothermia/Hypothermia-EEG-Data/'; 

sm=0;  ss=0.03;
CaseNo=1; 
for CaseNo=[1 2 3 4 5 6 7 8 10 12 13 ]; % [4 5 6 7 8 
    full_path = [td 'cs' num2str(CaseNo) '_bursts'];
    load(full_path)
%     str=['load cs' num2str(CaseNo) '_bursts']; eval(str); %  Burst Tburst T s z Fs
    
    if CaseNo==2; s=s(1216600:end); T=T(1216600:end); end
    if CaseNo==4; s=s(1211000:4675402); T=T(1211000:4675402); end
    if CaseNo==5; s=s(28000:end); T=T(28000:end); end
    if CaseNo==6; s=s(840000:3568802); T=T(840000:3568802); end; %3568402); end
    if CaseNo==7; s=s(294000:end); T=T(294000:end); end
    if CaseNo==8; s=s(681800:end); T=T(681800:end); end
    if CaseNo==10; 
        load([td 'case8ind.mat'])
%         load case8ind;
        s(ind)=[]; T(ind)=[]; 
    end
    Nt=length(s); Fs=200; dt=1/Fs; t=(0:Nt-1)*dt; 
    movingwin=[.25 .01]; 
    s=locdetrend(s',Fs,movingwin)'; 
    s=s/std(s); 
    % Find envelope: (a) HPF, (b) square, (c) LPF
    tt=-5:dt:5; 
    sig=.3; g=exp(-1/2/sig^2*tt.^2); g=g/sum(g); fs=conv(s,g,'same'); fb=s-fs; fs=(fb).^2; 
    % Low pass the squared HP component
    sig=.08; ge=exp(-1/2/sig^2*tt.^2); ge=ge/sum(ge); 
    e=sqrt(conv(fs,ge,'same'));         

    %% Find stretches above thresh at least 0.5 sec; off time: <thresh at least 0.5 sec
    
    z=(e<.15);
    if CaseNo==5; z=(e<.2); end
    if CaseNo==6; z=(e<1.5); end
    if CaseNo==8; z=(e<.2); end
    figure(1); clf; 
    subplot(411); plot(t,s); axis tight
    subplot(412); imagesc(z); colormap gray
    
    zz=z(1:2:end); 
    binomial_trials = 200*10; tstart = tic; Fs=200; 
    [Responses, ptile,xnew,signewsq,A,p, x, ss] = BSP_EM_estim(zz,Fs,binomial_trials);
    tz=linspace(min(t),max(t),length(ptile.pmode)); 
    lower=ptile.p050; upper=ptile.p950;
    bsp=ptile.pmode; 
    
    subplot(413); plot(tz,bsp,'k'); 
    y=T; x=t; xi=tz; yi=interp1(x,y,xi); 
    subplot(414); plot(xi,yi,'k'); 
    
%     str=['save PlotDataForCaseNo' num2str(CaseNo) ' t s xi yi tz bsp z']; eval(str); 
    full_path = [td 'PlotDataForCaseNo' num2str(CaseNo)];
    save(full_path, 't', 's', 'xi', 'yi', 'tz', 'bsp', 'z')
    
end



