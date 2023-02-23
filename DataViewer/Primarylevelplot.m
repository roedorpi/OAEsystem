function [handles] = Primarylevelplot(handles)
ax2 = handles.SubFigL;
ax = [800 5500 20 100];
set(ax2, 'Visible','on')
set(get(ax2,'XLabel'),'String','Frequency [Hz]')
set(get(ax2,'YLabel'),'String','dB SPL re. 20\muPa')
set(ax2,'Xlim',ax(1:2),'ylim',ax(3:4),'Box','on',...
    'YGrid','on','XGrid','on');
if ~handles.holdplot
    cla(ax2)
end
set(ax2,'Nextplot','add');
handles.plplot(1) = plot(ax2,handles.Freqs(:,1),handles.Pri_lev(:,1),...
    char(handles.coll(handles.numberofplots)),'linewidth',2, ...
    'markersize',6,'markerfacecolor',...
    char(handles.colll(handles.numberofplots)));
handles.plplot(2) = errorbar(ax2,handles.Freqs(:,1),handles.Pri_lev(:,1), ...
    handles.Pri_lev(:,2),char(handles.col(handles.numberofplots)));
handles.plplot(3) = plot(ax2,handles.Freqs(:,2),handles.Pri_lev(:,3),...
    char(handles.coll(handles.numberofplots)),'linewidth',2, ...
    'markersize',6,'markerfacecolor',...
    char(handles.colll(handles.numberofplots)));
handles.plplot(4) = errorbar(ax2,handles.Freqs(:,2),handles.Pri_lev(:,3), ...
    handles.Pri_lev(:,4),char(handles.col(handles.numberofplots)));
text(ax(1), 65, 'L1','Fontsize',14,'parent',ax2)
text(ax(1), 45, 'L2','Fontsize',14,'parent',ax2)
% freqs = 0:1/handles.T:handles.fs-1;
%     for k = 1:16
%         set(ax2,'Nextplot','add');
%         plot(ax2,freqs(1:length(TS)/2),20*log10(TS_f(1:length(TS)/2,k:16:end)./20e-6))
%         plot(ax2,F1(k,1),20*log10(f1_lev(k,:)./20e-6),'*','markersize',8)
%         plot(ax2,F2(k,1),20*log10(f2_lev(k,:)./20e-6),'s','markersize',8)
%         plot(ax2,FDP(k,1),20*log10(dplev(k,:)./20e-6),'o','markersize',8)
%       
%         M(k) = getframe(ax2);
%         cla(ax2);
% 
%     end
%     movie(ax2,M,10,3);
%    guidata(hObject, handles);

end