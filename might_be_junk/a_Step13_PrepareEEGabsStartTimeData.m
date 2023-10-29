clear all; clc; format compact

%% Load EEG data. Temp comes with it

% td='..\..\BurstSuppressionPaper_Material\ForMP\Hypothermia';
td = '..//Hypothermia-EEG-Data';

pd=pwd; 

fileNo=[1 2 3 4 5 7 8 10 12 13 14];
for i=fileNo
    disp(i); 
    cd(td); eval(['load h' num2str(i) '_data']); cd(pd); 
    EEGSTARTTIME_ABS(i)=ts(1); 
end

save EEGSTARTTIME EEGSTARTTIME_ABS