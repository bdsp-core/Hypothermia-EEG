function [t,s,en,Fs,EEGSTARTTIME_ABS,th,smallTh,uevents,ueventtimes,t_temp,Temp]=fcnGetDetrendedEEG(fileNo);

% tsc='..\..\BurstSuppressionPaper_Material\SpectralSegmentation\Score BSup';
tsc='../Hypothermia-EEG-Data/'
% td='..\..\BurstSuppressionPaper_Material\SpectralSegmentation\Hypothermia';
td = tsc; 
pd=pwd; 

%% Load EEG data. Temp comes with it
% cd(td); 
full_path = [tsc 'h' num2str(fileNo) '_data'];
load(full_path);
% eval(['load h' num2str(fileNo) '_data']); 
% cd(pd); 
EEGSTARTTIME_ABS=ts(1); 
t=(1:length(data(1,:)))/Fs; % time relative to beginning of EEG recording, in seconds
data=-data(1,:);  data=data-mean(data); dt=1/Fs; Nt=length(data); s=data; 

t_temp=(t_temp-EEGSTARTTIME_ABS)*24*60*60; % rezero and convert to seconds
Temp=T; 

%% get rid of high freq noise
tt=linspace(-1,1,151); sig=.1; g=exp(-(tt/sig).^2); g=g/sum(g); s=conv(s,g,'same'); dd=s; 
%% detrend 
dt= (t(2)-t(1)); Fs=1/dt; Nt=10*Fs; sig=2; tt=linspace(-5,5,Nt); g=exp(-(tt/sig).^2); g=g/sum(g); 
lp=conv(s,g,'same'); s=s-lp;

%% Get envelope
nt=Fs*3; ttt=linspace(-1.5,1.5,nt); sig=.6; g=exp(-(ttt/sig).^2); g=g/sum(g); 
en=sqrt(conv(s.^2,g,'same')); 

%% For each EEG need a threshold (mean + 3SD of background)
cd(tsc); eval(['load h' num2str(fileNo) '_ScoringLM']); cd(pd); % uevents ueventtimes stages stagetimes
[~,jj]=sort(ueventtimes(:,1)); uevents=uevents(jj); ueventtimes=ueventtimes(jj,:); % Put events in order -- sometimes they are not

ind=[];
for i=1:length(uevents); 
    if uevents(i)==2; 
        ind=[ind find(t>=ueventtimes(i,1) & t<=ueventtimes(i,2))]; 
    end
end
m=mean(en(ind)); sig=std(en(ind)); th=m+2*sig; smallTh=m; 
