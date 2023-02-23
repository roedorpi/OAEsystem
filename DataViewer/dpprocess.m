function [handles] = dpprocess(filename, handles)
%read from filename
ID = filename(1:3);
Ear = filename(8);
MeasTime =  datenum(regexp(filename,...
    '\d\d-\d\d-\d\d\d\d\ \d\d\-\d\d-\d\d','match'),'dd-mm-yyyy HH-MM-SS');
handles.MainFigLegend = strcat(ID,'-', Ear,'@',datestr(MeasTime));
%prepare axes
ax1 = handles.MainFig;
set(ax1, 'Visible','on','Fontsize',14)
set(get(ax1,'XLabel'),'String','Frequency [Hz]','Fontsize',14)
ax = [700 3500 -30 25];
set(ax1,'Xlim',ax(1:2),'ylim',ax(3:4),...
    'Box','on','XGrid','on','YGrid','on');
k = handles.numberofplots;
n = handles.holdplot;
if ~k && ~n
    k = 1;
    cla(ax1)
elseif ~k && n
    k = 1;
elseif k>0 && ~n 
    k = 1;
    cla(ax1)
elseif k>0 && n
    k = k+1;
end

if k <= length(handles.col)
    handles.numberofplots = k;
    set(ax1,'nextplot','add');
    plot(ax1,handles.Freqs(:,3),handles.Dp_lev(:,1),char(handles.coll(k)),...
        'linewidth',2,'markersize',8,'markerfacecolor', char(handles.colll(k)))
    plot(ax1,[2500 2600],[24-(k*3) 24-(k*3)],char(handles.col(k)),'linewidth',2)
    text(2610,24-(k*3), handles.MainFigLegend,'parent',ax1)
    text(ax(1), ax(4)+2,'dB SPL re. 20\muPa','fontsize',14,'parent',ax1)
    
    
    if handles.DPstdplot
        handles = DPstdplot(handles);   
    end
    if handles.noisefloorplot
        handles = Noisefloorplot(handles);
    end
    if handles.Primarylevelplot
        handles = Primarylevelplot(handles);
    end
    if handles.Primaryratioplot
        handles = Primaryratioplot(handles);
    end

else
    errordlg('Too many plots in the graph!! Clear the axes.','Plotting error','modal');
end



end













%[H,P,CI,STATS] = ttest(DPLevel(:,:,1),DPLevel(:,:,end),0.05)


% clearvars -except dur fs num_points num_reps timesignals
% 
% freqs = 0:1/dur:fs/2;
% freqdomain = 20*log10(2*abs(fft(timesignals'))/(dur*fs)/20e-6);
% freqdomain = freqdomain(1:dur*fs/2+1,:);
% size(freqdomain)
% semilogx(freqs,freqdomain(:,2))
% xlim([1 fs/2])
% grid on





