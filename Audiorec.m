%% Example Session 5 Introduction to MATLAB
% Audio recorder
%
% Author: Rodrigo Ordoñez <rop@es.aau.dk>
%%
%% Create GUI Figure:

function Audiorec
%%
hmainfig = figure(...
    'MenuBar','none',...
    'name', 'AUDIO RECORDER',...
    'Position', [200 300 1024 500],...
    'resize', 'on',...
    'visible', 'on',...
    'tag', 'mainfig'); % 'toolbar', 'none',...

ax1 = axes('Parent',hmainfig,'position',[0.07 0.6 0.6 0.35],'tag','Ax1');
xlabel(ax1, 'Time [s]');
ylabel(ax1, 'Amplitude');

ax2 = axes('Parent',hmainfig,'position',[0.07 0.12 0.6 0.35],'tag','Ax2');
xlabel(ax2, 'Frequency [Hz]');
ylabel(ax2, 'Amplitude [dB]');

hsettings = uipanel('Parent',hmainfig,'Title','Controls',...
             'Position',[.72 .35 .25 .6]);

hdysplay = uipanel('Parent',hmainfig,'Title','Audio info',...
             'Position',[.72 .1 .25 .2]);
         
hstart = uicontrol(hsettings,'style','togglebutton',...
    'string','Start/Stop', 'Value', 0,'units','normalized',...
    'position', [0.1 0.85 0.8 0.15],'CallBack',{@RecordAudio},...
    'tag','Start');

hperiod = uicontrol(hsettings,'style','slider',...
    'units','normalized','position', [0.1 0.65 0.8 0.05],...
    'Min',0.02,'Max',5,'SliderStep',[0.001 0.05],'value',5,...
    'CallBack',{@TimePeriod},'Tag','Period');

hst1 = uicontrol(hsettings,'Style','text',...
    'String',sprintf('Analysis time: %1.3f seconds',get(hperiod,'value')),...
    'units','normalized','Position',[0 0.7 1 0.1],'tag','st1'); 

hchan = uicontrol(hsettings,'Style','popup','string','1|2|[1 2]',...
    'value',1,'callback',{@ChanSelect},...
    'units','normalized','Position',[0.7 0.5 0.2 0.1],'tag','chan');

hst4 = uicontrol(hsettings,'Style','text',...
    'String','Recording Channel',...
    'units','normalized','Position',[0 0.5 0.6 0.1],'tag','st4'); 

hcalibrate = uicontrol(hsettings,'style','pushbutton',...
    'string','Calibrate','units','normalized',...
    'position', [0.1 0.38 0.8 0.1],'CallBack',{@Calibrate},...
    'tag','Calib');

hcaldB = uicontrol(hsettings,'Style','popup',...
    'string','94dB@1k|95.2dB@1k|123.3dB@250|123.8dB@250',...
    'value',1,'callback',{@CaldBSelect},...
    'units','normalized','Position',[0.3 0.27 0.4 0.1],'tag','chan');


hst2 = uicontrol(hsettings,'Style','text',...
    'String','Not calibrated',...
    'units','normalized','Position',[0 0.01 1 0.2],'tag','st2'); 

hst3 = uicontrol(hdysplay,'Style','text',...
    'String','',...
    'units','normalized','Position',[0 0.1 1 0.8],'tag','st3'); 


Data = struct(...
    'FS', 48000,...
    'T', get(hperiod,'value'),...
    'AudioData', [],...
    'Correction', 1,...
    'Chan',get(hchan,'Value'),...
    'CaldBSelect', 94);


handles = guihandles(hmainfig);
handles.Data = Data;
guidata(hmainfig,handles);
end

function TimePeriod(Obj,Event)
    handles = guidata(Obj);
    handles.Data.T = get(Obj,'value');
    set(handles.st1,'String',...
        sprintf('Analysis time: %1.3f seconds',get(Obj,'value')));
    guidata(Obj,handles)
end
    
function SampleFrequency(Obj,Event)
    handles = guidata(Obj);
    str = get(Obj,'String');
    val = get(Obj,'Value');
    handles.Data.FS = str2num(str{val});
    if ~isempty(handles.Data.AudioRec)
       stop(handles.Data.AudioRec)
       handles = RecInit(handles);
       set(handles.StartStop,'value',0)
    end
    guidata(Obj,handles)
end

function ChanSelect(Obj,Event)
    handles = guidata(Obj);
    ChanSel = get(Obj,'Value');
    switch ChanSel
        case 1
            Chan = 1;
        case 2 
            Chan = 2;
        case 3
            Chan = [1 2];
    end
    handles.Data.Chan = Chan; 
    guidata(Obj,handles)
    
end

function CaldBSelect(Obj,Event)
    handles = guidata(Obj);
    RefVal = get(Obj,'Value');
    switch RefVal
        case 1
            dBVal = 94;
        case 2
            dBVal = 95.2;
        case 3
            dBVal = 123.3;    
        case 4
            dBVal = 123.8;
    end
    handles.Data.CaldBSelect = dBVal;
    guidata(Obj,handles)
end

function RecordAudio(Obj,Event)
    handles = guidata(Obj);
    ButtonPos = get(Obj,'value');
    if ButtonPos 
        AudioData = handles.Data.AudioData;
        AutoPlayrecInit
        while ButtonPos
            set(handles.Period,'enable','off')
            RecLength = floor(handles.Data.T*handles.Data.FS);

            out = playrec('rec',RecLength,handles.Data.Chan);
            playrec('block',out);
            Out = playrec('getRec',out);

            AudioData = [AudioData; Out]; 
            PlotAudio(AudioData,handles);
            ButtonPos = get(Obj,'value');
        end
        handles.Data.AudioData = AudioData;
    else
        playrec('reset')
        set(handles.Period,'enable','on')
    end
    
    guidata(Obj,handles)
end


function Calibrate(Obj,Event)
    handles = guidata(Obj);
    AutoPlayrecInit
    
    RecLength = floor(2*handles.Data.FS);
    out = playrec('rec',RecLength,handles.Data.Chan);
    playrec('block',out);
    y = playrec('getRec',out);
    playrec('reset')
    
    RMS = sqrt(sum(y.^2)/length(y));
    y_dB = 10*log10(RMS.^2/4e-10);
    Correction = 10^((handles.Data.CaldBSelect - y_dB)/20);
    y = y.*Correction; %Time signal in Pa. 

    RMS_ = sqrt(y'*y/length(y));
    y_dB_ = 10*log10(RMS_^2/4e-10);
    
    CalInfo = sprintf('Correction Factor: %4.5f \n Measured SPL = %4.1f dB \n Adjusted SPL %4.1f', Correction, y_dB, y_dB_);
    set(handles.st2,'string',CalInfo);
    handles.Data.Correction = Correction;
    guidata(Obj,handles)
end

%% PlotAudio: TimerFcn CallBack
% Mandatory input arguments set from the calling function: 
% 
% obj: handle to in/out object
%
% Event: When should the calling function execute this script
%
% Optional arguments can be set from the calling function and are available 
% through the varargin cell array.
%
function PlotAudio(AudioData,handles)
%% 
% get audiorecorder information.
FS = handles.Data.FS;
T = handles.Data.T;

TotSamples = length(AudioData);
CurrentSamples = AudioData(end-floor(FS*T)+1:end,:);

%%
% Get the axes handles set in the |UserData| field of the |audiorecorder|
% defined in the |CreateFigure| CallBack.

ax1 = handles.Ax1;
ax2 = handles.Ax2;
%%
% calculate number of samples in each callback period and make a time
% vector
NrSamples = floor(T*FS);
t = 1/FS:1/FS:T;
%%
% Get the audio data, as double between -1 and 1. Keep only the last T 
% seconds of the recording, and calculate the RMS value.         
y = CurrentSamples;
% convert to pascals
Correction = handles.Data.Correction;
y = y.*Correction; %Time signal in Pa. 
% Calculate RMS and dB.
RMS_ = sqrt(sum(y.^2)/NrSamples); 
y_dB_ = 10*log10(RMS_.^2/4e-10);

%%
% plot in the time domain
plot(ax1,t,y,[t(1) t(end)],[RMS_' RMS_'],'linewidth',2);
set(ax1, 'xlim',[0 T],'Ylim',[-3 3])

%%
% Make the frequency analysis with the fft fuction (Fast Fourier Transform)
% and create a frequency vector for plotting.
% 10*log10((2.*abs(fft(signal))./length(signal)).^2./4e-10);
Y = 10*log10((2.*abs(fft(y))./length(y)).^2./4e-10);
F = 0:FS/NrSamples:FS-FS/NrSamples;
%%
% plot the frequency amplitude in dB re. 20 muPa.  

semilogx(ax2, F(1:round(end/2)),Y(1:round(end/2),:))
set(ax2, 'xlim',[20 FS/2], 'ylim',[0,110])

%%
% Write on the figure how many seconds have been recorded:
RecSecs = TotSamples/FS;
AudioInfo = sprintf('Total Running Time %10.2f seconds\n SPL = %4.1f dB re. \n', RecSecs, y_dB_);
set(handles.st3,'string',AudioInfo);
drawnow
end
