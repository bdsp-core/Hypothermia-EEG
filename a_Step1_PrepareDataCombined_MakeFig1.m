clear all; clc; format compact;  

%% Creates event files, binary vector files, spectra files
td = '../Hypothermia-EEG-Data/';
% These specify the useful part of the data
% Note time convention: refer all to eeg elapsed time;
fileNo=[1 2 3 4 5 7 8 10 12 13 14];
Tstart(fileNo)=[0   105 0   105 5   35  60  20  60  0   60]*60; % in sec
Tend(fileNo)=  [500 300 260 500 150 250 500 180 190 200 140]*60; % in sec

for fileNo=[1 2 3 4 5 7 8 10 12 13 14] % [1 2 3 4 5 7 8 10 12 13 14]
    disp(fileNo);
    
    %% Creates event files, binary vector files, spectra files [h_ScoringXX, BinaryXX, EVENTINDSXX, SPECTRAXX]
    %% Automatic segmentation of hypothermia files 
    %*****************************************
    %% Creates event files, binary vector files, spectra files
    pd=pwd;
    
    %% Get detrended EEG and temperature
    [t,s,en,Fs,EEGSTARTTIME_ABS,th,smallTh,uevents,ueventtimes,t_temp,Temp]=fcnGetDetrendedEEG(fileNo);     %% Find onset and offset times for bursts, save event file
    %% Get event data
    [uevents,ueventtimes,ei,ARTINDS,BURSTINDS,BURST]=fcn_SegmentEEG(Fs,t,en,th,uevents,ueventtimes); 
    %% Extract bursts and suppressions for analyses
    [SPECTRA,BurstWaveform,tb,bl,bt,f]=fcnGetIsolatedBursts(t,s,Fs,en,th,t_temp,Temp,smallTh,BURST);
    
    %% throw out bad data at ends
    ind=find(t<=Tstart(fileNo) | t>Tend(fileNo)); tz=t; tz(ind)=[]; ei(ind)=[];
    
    %% Get binary vector for BSR calculation
   [tbsp,bsp,z]=fcnGetBSP(tz,ei,fileNo);
    
    %% Save data
%     tsc='..\..\BurstSuppressionPaper_Material\SpectralSegmentation\Score BSup'; pd=pwd;
     stages=[]; stagetimes=[]; 
%     cd(tsc); 
%     eval(['save h' num2str(fileNo) '_ScoringAuto uevents ueventtimes stages stagetimes']); cd(pd); 
%     eval(['save CURATEDDATA' num2str(fileNo) ' t Fs s z tz ei en th t_temp Temp tbsp bsp ARTINDS BURSTINDS BURST']); 
%     eval(['save SPECTRA' num2str(fileNo) ' SPECTRA BurstWaveform tb bl bt f']);
%     eval(['save INDIVBURSTDATA' num2str(fileNo) ' t bt s Fs en th smallTh BURST']); 

    eval(['save ''' td 'h' num2str(fileNo) '_ScoringAuto''' 'uevents ueventtimes stages stagetimes']);
    eval(['save ''' td 'CURATEDDATA' num2str(fileNo) ''' t Fs s z tz ei en th t_temp Temp tbsp bsp ARTINDS BURSTINDS BURST']);
    eval(['save ''' td 'SPECTRA' num2str(fileNo) ''' SPECTRA BurstWaveform tb bl bt f']);
    eval(['save ''' td 'INDIVBURSTDATA' num2str(fileNo) ''' t bt s Fs en th smallTh BURST']);  
   
    %% Make figure
    fcnPlotMainFig(fileNo)
end

