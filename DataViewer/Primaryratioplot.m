function [handles] = Primaryratioplot(handles)
ax3 = handles.SubFigR;
set(ax3, 'Visible','on')
set(get(ax3,'XLabel'),'String','Frequency [Hz]')
set(get(ax3,'YLabel'),'String','f_2/f_1')
set(ax3,'Xlim',[800 3500])
set(ax3, 'NextPlot','replacechildren','YAxisLocation','right')
handles.Priratioplot = plot(ax3,handles.Freqs(:,3),...
    handles.Freqs(:,2)./handles.Freqs(:,1),...
    char(handles.col(handles.numberofplots)),'linewidth',2);
set(ax3,'Box','on' ,'YGrid','on','XGrid','on');
%    guidata(hObject, handles);
end