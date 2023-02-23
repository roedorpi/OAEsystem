function [handles] = Noisefloorplot(handles)
%plot noise floor of DPOAE
Noise = cat(2, handles.Dp_lev(:,3)-handles.Dp_lev(:,4),...
    handles.Dp_lev(:,3),handles.Dp_lev(:,3)+handles.Dp_lev(:,4));
Y = cat(1,Noise(:,3),Noise(end:-1:1,1));
%Y = cat(1,Noise(:,1),Noise(end:-1:1,2));
X = cat(1,handles.Freqs(:,3),handles.Freqs(end:-1:1,3));
Z = ones(2*length(handles.Freqs(:,3)),1);

Gray = [0.6 0.6 0.6
        0.7 0.7 0.7
        0.8 0.8 0.8
        0.9 0.9 0.9];
    
Pos = get(handles.MainFig, 'Position');
Fig = get(handles.MainFig, 'Parent');
Hnp = axes('Parent',Fig,'Position',Pos);
if handles.numberofplots < 4
handles.p1(handles.numberofplots) = patch(handles.MainFig,X,Y,'b');
set(handles.p1(handles.numberofplots),...
    'facecolor',Gray(handles.numberofplots,:),'Edgecolor','none','facealpha',0.5);
set(Hnp,'box','on');
else
    handles.p1(handles.numberofplots) = 1;
end
handles.p2(1,handles.numberofplots) = ...
    plot(handles.MainFig,handles.Freqs(:,3),handles.Dp_lev(:,3),...
    char(handles.colll(handles.numberofplots)),'linewidth',2);
handles.p2(2,handles.numberofplots) = ...
    plot(handles.MainFig,handles.Freqs(:,3),handles.Dp_lev(:,3) ...
    + handles.Dp_lev(:,4),char(handles.colll(handles.numberofplots)));
handles.p2(3,handles.numberofplots) = ...
    plot(handles.MainFig,handles.Freqs(:,3),handles.Dp_lev(:,3) ...
    - handles.Dp_lev(:,4),char(handles.colll(handles.numberofplots)));
%guidata(hObject, handles); 

end

% Noise = cat(2, varargin{2}.Dp_lev(:,3)-varargin{2}.Dp_lev(:,4),...
%     varargin{2}.Dp_lev(:,3),varargin{2}.Dp_lev(:,3)+varargin{2}.Dp_lev(:,4));
% Y = cat(1,Noise(:,3),Noise(end:-1:1,1));
% X = cat(1,varargin{2}.Freqs(:,3),varargin{2}.Freqs(end:-1:1,3));
% Z = ones(2*length(varargin{2}.Freqs(:,3)),1);
% fill(X,Y,Z,'facecolor',[0.8 0.8 0.8],'Edgecolor','none');
% plot(ax1,varargin{2}.Freqs(:,3),varargin{2}.Dp_lev(:,3),'k-')
% plot(ax1,varargin{2}.Freqs(:,3),varargin{2}.Dp_lev(:,3) + varargin{2}.Dp_lev(:,4),'k-')
% plot(ax1,varargin{2}.Freqs(:,3),varargin{2}.Dp_lev(:,3) - varargin{2}.Dp_lev(:,4),'k-')
