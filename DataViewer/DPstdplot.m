function [handles] = DPstdplot(handles) 
handles.wiskers(handles.numberofplots) = errorbar(handles.MainFig,...
    handles.Freqs(:,3),handles.Dp_lev(:,1), ...
    handles.Dp_lev(:,2),char(handles.col(handles.numberofplots)));
%guidata(hObject, handles); data is updates through the return variable
end