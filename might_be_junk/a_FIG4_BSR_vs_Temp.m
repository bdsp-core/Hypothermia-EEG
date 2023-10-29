clear all; clc; format compact; 

%% Compute spectra of suppressions in cardiac surgery hypothermia data
Fs=200; 
ct=0; PP=[]; PP0=[]; PP1=[]; BT=[]; % group ensemble    
td = '/Users/mbw/cdac Dropbox/Datasets/zz_Hypothermia/Hypothermia-EEG-Data/'; 

for CaseNo=[1 2 3 4 5 6 8 10 12 13  ]; % [4 5 6 7 8 10 12 13]
    disp(CaseNo); 
    %str=['load EEG_Z_TEMP_' num2str(CaseNo)]; eval(str); % get: ts s z T']; 

    filename = ['PlotDataForCaseNo' num2str(CaseNo)];
    full_path = [td filename];
    load(full_path);  % t s xi yi tz bsp z
    ts=t-t(1); 
    ts=ts/60/60; xi=xi/60/60;
    yi=(yi-32)/1.8;
    %% Fix times -- change to hours
    
    %% Figures
    L=.1; W=.8; B=.7; H=.25;
    figure(1); clf; 
    h=subplot('position',[L B W H]); set(h,'xtick',[]);  
    h=plot(ts,s,'k'); axis([0 max(ts) -20 20]); 
    box off; 
    set(gca,'xtick',[]); 
    set(gca,'FontName','Arial','FontSize',14); 
    f=figure(1); 
    set(gca,'xcolor','w','tickdir','out'); 
    ylabel('\muV','FontName','Arial','fontsize',14); 

    % Binary data
    B=.66; H=.02;
    h=subplot('position',[L B W H]); set(h,'xtick',[]); if ct>5; set(h,'ytick',[]);   end
    imagesc(z); colormap gray; axis off; 

    % BSP
    B=.38; H=.25;
    h=subplot('position',[L B W H]); set(h,'xtick',[]); 
    plot(xi,bsp,'k','linewidth',2); 
    ylabel('BSP','FontName','Arial','fontsize',14); 
    set(gca,'FontName','Arial','FontSize',14); 
    axis([0 max(xi) -.01 1.01]);  
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
    
    
     % Temp
    B=.1; H=.25;
    h=subplot('position',[L B W H]); 
    plot(xi,yi,'k','linewidth',2); 
    ylabel('Temperature (\circC)','FontName','Arial','fontsize',14); 
    xlabel('Time (hours)','FontName','Arial','fontsize',14); 

    set(gca,'FontName','Arial','FontSize',14); 
    axis([0 max(ts) 16 36]);  
    
    set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YTick'       , 0:.25:1,...
    'LineWidth'   , 1         );
    set(gcf,'color',[1 1 1]); 

    grid on; %grid minor
    set(gca,'ytick',[16:5:34])
    set(gca,'xtick',[0:0.5:10]); 
    %% Save figure
    figure(1);
    
    FileName=['Fig4_EEG_BSR_Temp_' num2str(CaseNo)]; 
    set(gcf, 'PaperPositionMode', 'auto'); 
    print('-dpng','-r300',FileName);
    
    %g=input('ok'); 
end
