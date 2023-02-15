%% Input and Output audio recorder using playrec
%
% Author: Rodrigo Ordoñez <rop@es.aau.dk>
%%
%% Create GUI Figure:

function AudioRecPlay
%%
addpath /home/rop/Documents/MATLAB/mls/
ScrN = get(0,'ScreenSize');

hmainfig = figure(...
    'MenuBar','none',...
    'name', 'AUDIO RECORDER',...
    'Position', ScrN,...   % [200 300 1024 500],...
    'resize', 'on',...
    'visible', 'on',...
    'tag', 'mainfig'); % 'toolbar', 'none',...

ax1 = axes('Parent',hmainfig,'position',[0.05 0.6 0.45 0.35],'tag','Ax1',...
    'box','on');
xlabel(ax1, 'Time [s]');
ylabel(ax1, 'Amplitude');
l1ax1 = line([0 1], [0 0],'color','b','parent',ax1,'visible','off');
l2ax1 = line([0 1], [0 0],'color','g','parent',ax1,'visible','off');
l3ax1 = line([0 1], [0 0],'color','r','parent',ax1,'visible','off');
l4ax1 = line([0 1], [0 0],'color','c','parent',ax1,'visible','off');



ax2 = axes('Parent',hmainfig,'position',[0.05 0.14 0.45 0.35],'tag','Ax2',...
    'box','on','xscale','log');
xlabel(ax2, 'Frequency [Hz]');
ylabel(ax2, 'Amplitude [dB]');

l1ax2 = line([0 1], [0 0],'color','b','parent',ax2,'visible','off');
l2ax2 = line([0 1], [0 0],'color','g','parent',ax2,'visible','off');
l3ax2 = line([0 1], [0 0],'color','r','parent',ax2,'visible','off');
l4ax2 = line([0 1], [0 0],'color','c','parent',ax2,'visible','off');


% control panels
setin = uipanel('Parent',hmainfig,'Title','Input Settings',...
             'Position',[.52 .35 .14 .6]);
setout = uipanel('Parent',hmainfig,'Title','Output Settings',...
             'Position',[.68 .35 .14 .6]);
setmeas = uipanel('Parent',hmainfig,'Title','Measurement',...
             'Position',[.84 .22 .14 .73]);
hdysplay = uipanel('Parent',hmainfig,'Title','Audio info',...
             'Position',[.52 .1 .3 .23]);
% Start Button
hstart = uicontrol(hmainfig,'style','togglebutton',...
    'string','Start/Stop', 'Value', 0,'units','normalized',...
    'position', [0.84 0.1 0.14 0.1],'CallBack',{@RecordAudio},...
    'tag','Start');

% Measurement set-up 
hperiod = uicontrol(setmeas,'style','slider',...
    'units','normalized','position', [0.1 0.85 0.8 0.05],...
    'Min',0.02,'Max',5,'SliderStep',[0.001 0.05],'value',mean([0.02 5]),...
    'CallBack',{@TimePeriod},'Tag','Period');

hst1 = uicontrol(setmeas,'Style','text',...
    'String',sprintf('Analysis time: %1.3f s',get(hperiod,'value')),...
    'units','normalized','Position',[0 0.9 1 0.05],'tag','st1'); 

uicontrol(setmeas,'Style','text',...
    'String','Stimulus type:',...
    'units','normalized','Position',[0 0.7 1 0.05]); 

stim = uicontrol(setmeas,'Style','popup','string','tone|impulse|MLS',...
    'value',1,'callback',{@StimSelect},...
    'units','normalized','Position',[0.1 0.65 0.8 0.05],'tag','stim');

% MLS measurements
uicontrol(setmeas,'Style','text', 'visible','off',...
    'String','Pre-averages:',...
    'units','normalized','Position',[0 0.55 0.6 0.05],'tag','nravertext'); 
hnraver = uicontrol(setmeas,'style','edit',...
    'units','normalized','position', [0.61 0.55 0.3 0.05],...
    'string',8,'Tag','nraver', 'visible','off');
uicontrol(setmeas,'Style','text', 'visible','off',...
    'String','Repetitions:',...
    'units','normalized','Position',[0 0.45 0.6 0.05],'tag','nrrepstext'); 
hnrreps= uicontrol(setmeas,'style','edit', 'visible','off',...
    'units','normalized','position', [0.61 0.45 0.3 0.05],...
    'string',4,'Tag','nrreps');
uicontrol(setmeas,'Style','text', 'visible','off',...
    'String','MLS order:',...
    'units','normalized','Position',[0 0.35 0.6 0.05],'tag','mlsordertext'); 
horder = uicontrol(setmeas,'style','edit', 'visible','off',...
    'units','normalized','position', [0.61 0.35 0.3 0.05],...
    'string',12,'Tag','mlsorder');

% Tone set up
uicontrol(setmeas,'Style','text', 'visible','on',...
    'String','Tone duration:',...
    'units','normalized','Position',[0 0.55 0.6 0.05],'tag','toneperiodtext'); 
htoneperiod = uicontrol(setmeas,'style','edit',...
    'units','normalized','position', [0.61 0.55 0.3 0.05],...
    'string',0.02,'Tag','toneperiod', 'visible','on');

uicontrol(setmeas,'Style','text', 'visible','on',...
    'String','EQ filter:',...
    'units','normalized','Position',[0 0.45 0.6 0.05],'tag','EQfilttext'); 
hnrreps= uicontrol(setmeas,'style','checkbox', 'visible','on',...
    'units','normalized','position', [0.61 0.45 0.3 0.05],...
    'value',0,'Tag','EQfilt');


% Input set-up
% Input 1
uicontrol(setin,'Style','text',...
    'String','Input 1:',...
    'units','normalized','Position',[0 0.9 0.6 0.05]); 

hch1 = uicontrol(setin,'Style','popup','string','on|off',...
    'value',1,'callback',{@ChanSelect},...
    'units','normalized','Position',[0.61 0.9 0.2 0.05],'tag','chan1');

% Input 2
uicontrol(setin,'Style','text',...
    'String','Input 2:',...
    'units','normalized','Position',[0 0.8 0.6 0.05]); 

hch2 = uicontrol(setin,'Style','popup','string','on|off',...
    'value',1,'callback',{@ChanSelect},...
    'units','normalized','Position',[0.61 0.8 0.2 0.05],'tag','chan2');

% Calibrate input
uicontrol(setin,'Style','text',...
    'String','Calib. Ref.',...
    'units','normalized','Position',[0 0.6 0.58 0.05]); 

hcaldB = uicontrol(setin,'Style','popup',...
    'string','94dB@1k|95.2dB@1k|123.3dB@250|123.8dB@250',...
    'value',1,'callback',{@CaldBSelect},...
    'units','normalized','Position',[0.59 0.6 0.4 0.05],'tag','cal');

hcalibrate = uicontrol(setin,'style','pushbutton',...
    'string','Calibrate','units','normalized',...
    'position', [0.1 0.38 0.8 0.1],'CallBack',{@Calibrate},...
    'tag','Calib');

uicontrol(setin,'Style','text',...
    'String','Input 1 not calibrated',...
    'units','normalized','Position',[0 0.15 1 0.13],'tag','caltext1'); 
uicontrol(setin,'Style','text',...
    'String','Input 2 not calibrated',...
    'units','normalized','Position',[0 0.01 1 0.13],'tag','caltext2'); 


% Audio info
uicontrol(hdysplay,'Style','text',...
    'String','Total Running Time: ',...
    'units','normalized','Position',[0 0.85 1 0.1],'tag','runtime'); 

uicontrol(hdysplay,'Style','text',...
    'String','Input 1:',...
    'units','normalized','Position',[0.01 0.01 0.48 0.7],'tag','info1'); 

uicontrol(hdysplay,'Style','text',...
    'String','Input 2:',...
    'units','normalized','Position',[0.51 0.01 0.48 0.7],'tag','info2'); 

% Output set-up
% Output 1
uicontrol(setout,'Style','text',...
    'String','Output 1:',...
    'units','normalized','Position',[0 0.9 0.6 0.05]); 

hout1 = uicontrol(setout,'Style','popup','string','on|off',...
    'value',2,'callback',{@ChanSelect},...
    'units','normalized','Position',[0.61 0.9 0.2 0.05],'tag','out1');

hlev1 = uicontrol(setout,'style','edit',...
    'units','normalized','position', [0.61 0.65 0.3 0.05],...
    'string',0.5,'CallBack',{@LevelAdjust},'Tag','Lev1');

uicontrol(setout,'Style','text',...
    'String','Digital Output 1:',...
    'units','normalized','Position',[0 0.65 0.6 0.05],'tag','Lev1text'); 

hfreq1 = uicontrol(setout,'style','edit',...
    'units','normalized','position', [0.61 0.55 0.3 0.05],...
    'string',1000,'CallBack',{@FreqAdjust},'Tag','Freq1');

uicontrol(setout,'Style','text',...
    'String','Frequency 1: ',...
    'units','normalized','Position',[0 0.55 0.6 0.05],'tag','Freq1text'); 

% Output 2
uicontrol(setout,'Style','text',...
    'String','Output 2:',...
    'units','normalized','Position',[0 0.8 0.6 0.05]); 

hout2 = uicontrol(setout,'Style','popup','string','on|off',...
    'value',2,'callback',{@ChanSelect},...
    'units','normalized','Position',[0.61 0.8 0.2 0.05],'tag','out2');

hlev2 = uicontrol(setout,'style','edit',...
    'units','normalized','position', [0.61 0.45 0.3 0.05],...
    'string',0.5,'CallBack',{@LevelAdjust},'Tag','Lev2');

uicontrol(setout,'Style','text',...
    'String','Digital Output 2:',...
    'units','normalized','Position',[0 0.45 0.6 0.05],'tag','Lev2text'); 

hfreq2 = uicontrol(setout,'style','edit',...
    'units','normalized','position', [0.61 0.35 0.3 0.05],...
    'string',1000,'CallBack',{@FreqAdjust},'Tag','Freq2');

uicontrol(setout,'Style','text',...
    'String','Frequency 2: ',...
    'units','normalized','Position',[0 0.35 0.6 0.05],'tag','Freq2text'); 

Data = struct(...
    'FS', 48000,...
    'T', get(hperiod,'value'),...
    'AudioData', [],...
    'Correction', [1 1],...
    'Chan1',1,...       % 1 in/out-put is on 
    'Chan2',1,...       % 0 in/out-put is off
    'Out1',0,...        % default both inputs are on and 
    'Out2',0,...        % both outputs are off
    'CaldBSelect', 94,...
    'OutLev1', str2double(get(hlev1,'string')),...
    'OutLev2', str2double(get(hlev2,'string')),...
    'OutFreq1', str2double(get(hfreq1,'string')),...
    'OutFreq2', str2double(get(hfreq2,'string')),...
    'NrAver', str2double(get(hnraver,'string')),...
    'NrReps', str2double(get(hnrreps,'string')),...
    'MLSorder', str2double(get(horder,'string')),...
    'TonePeriod', str2double(get(htoneperiod,'string')));
set(ax1, 'xlim',[0.02 Data.T],'Ylim',[-1 1])
set(ax2, 'xlim',[20 Data.FS/2],'Ylim',[-30 30])

handles = guihandles(hmainfig);
handles.Data = Data;

handles.l1ax1 = l1ax1;
handles.l2ax1 = l2ax1;
handles.l3ax1 = l3ax1;
handles.l4ax1 = l4ax1;
handles.l1ax2 = l1ax2;
handles.l2ax2 = l2ax2;
handles.l3ax2 = l3ax2;
handles.l4ax2 = l4ax2;

guidata(hmainfig,handles);
end

function TimePeriod(Obj,Event)
    handles = guidata(Obj);
    handles.Data.T = get(Obj,'value');
    set(handles.st1,'String',...
        sprintf('Analysis time: %1.3f seconds',get(Obj,'value')));
    guidata(Obj,handles)
end
function LevelAdjust(Obj,Event)
    handles = guidata(Obj);
    switch get(Obj,'Tag')
        case 'Lev1'
            handles.Data.OutLev1 = str2double(get(Obj,'string'));
        case 'Lev2'
            handles.Data.OutLev2 = str2double(get(Obj,'string'));
    end
    guidata(Obj,handles)
end    

function FreqAdjust(Obj,Event)
    handles = guidata(Obj);
    switch get(Obj,'Tag')
        case 'Freq1'
            handles.Data.OutFreq1 = str2double(get(Obj,'string'));
        case 'Freq2'
            handles.Data.OutFreq2 = str2double(get(Obj,'string'));
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
            Chan = 0;
    end
    switch get(Obj,'Tag')
        case 'chan1'
            handles.Data.Chan1 = Chan; 
        case 'chan2'
            handles.Data.Chan2 = Chan;
        case 'out1'
            handles.Data.Out1 = Chan;
        case 'out2'
            handles.Data.Out2 = Chan;
    end
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

function StimSelect(Obj,Event)
    handles = guidata(Obj);
    StimType = get(Obj,'value');
    StimLength = floor(handles.Data.T*handles.Data.FS);
    t = 0:1/handles.Data.FS:handles.Data.T-1/handles.Data.FS;
    Amplitude = [handles.Data.OutLev1 handles.Data.OutLev2];
    Freqs = [handles.Data.OutFreq1 handles.Data.OutFreq2];
    if handles.Data.Out1 && handles.Data.Out2
        PlayChan = [1 2];
    elseif handles.Data.Out1 && ~handles.Data.Out2
        PlayChan = 1;
    elseif ~handles.Data.Out1 && handles.Data.Out2
        PlayChan = 2;
    elseif ~handles.Data.Out1 && ~handles.Data.Out2
        PlayChan = [0 0];
    end
    switch StimType
        case 1
            set([handles.out1 handles.out2 handles.Lev1 handles.Lev2 ...
                handles.Freq1 handles.Freq2 handles.Period ...
                handles.chan1 handles.chan2], 'enable','on');
            
            set([handles.mlsorder handles.nraver handles.nrreps ...
                handles.mlsordertext handles.nravertext handles.nrrepstext],...
                'visible','off');
            set([handles.toneperiod handles.toneperiodtext...
                handles.EQfilt handles.EQfilttext],...
                'visible','on');

            set(handles.Start,'string','Start/Stop','CallBack',{@RecordAudio})
            T = str2double(get(handles.toneperiod,'string'));
            t = 0:1/handles.Data.FS:T-1/handles.Data.FS;
            winL=0.01; % 40ms length
            tt=0:1/handles.Data.FS:winL-1/handles.Data.FS;
            win=sin(2*pi*tt*(1/winL/2));
            window=[win(1:length(win)/2) ones(1,length(t)-length(win)) ...
                win((length(win)/2)+1:end)];
            if numel(PlayChan) == 2
                x(:,1) = Amplitude(1).*sin(2*pi*Freqs(1).*t).*window;
                x(:,2) = Amplitude(2).*sin(2*pi*Freqs(2).*t).*window;
                
            else
                if PlayChan == 1
                    x = Amplitude(1).*sin(2*pi*Freqs(1).*t).*window;
                else
                    x = Amplitude(2).*sin(2*pi*Freqs(2).*t).*window;
                end
                x = x';
            end
            if get(handles.EQfilt,'value');
                load EQFilt.mat
                if numel(PlayChan) == 2
                    x(:,1) = filter(hd1,x(:,1));
                    x(:,2) = filter(hd1,x(:,2));
                else
                    if PlayChan == 1
                        x = filter(hd1,x);
                    else
                        x = filter(hd2,x);
                    end
                
                end
            end
        case 2
            set([handles.out1 handles.out2 handles.Lev1 handles.Lev2 ...
                handles.Period handles.chan1 handles.chan2], 'enable','on');
            set([handles.mlsorder handles.nraver handles.nrreps ...
                handles.mlsordertext handles.nravertext handles.nrrepstext...
                handles.toneperiod handles.toneperiodtext],...
                'visible','off');
            set([handles.EQfilt handles.EQfilttext],...
                'visible','on');

            set(handles.Start,'string','Start/Stop','CallBack',{@RecordAudio})
            disp('Impulses!!!')
            
            Impulse = ones(1,round(80e-6*handles.Data.FS));
            Impulse = [zeros(1,100) Impulse ...
                zeros(1,round((20e-3-80e-6)*handles.Data.FS)-100)]; %impulse length 0.02 s
            
            if numel(PlayChan) == 2
                x(:,1) = Amplitude(1).*Impulse;
                x(:,2) = Amplitude(2).*Impulse;
            else
                if PlayChan == 1
                    x = Amplitude(1).*Impulse';
                else
                    x = Amplitude(2).*Impulse';
                end
            end
            
            if get(handles.EQfilt,'value');
                load EQFilt.mat
                if numel(PlayChan) == 2
                    x(:,1) = filter(hd1,x(:,1));
                    x(:,2) = filter(hd2,x(:,2));
                else
                    if PlayChan == 1
                        x = filter(hd1,x);
                    else
                        x = filter(hd2,x);
                    end
                
                end
            end
            
        case 3
            set([handles.out1 handles.out2 handles.Lev1 handles.Lev2 ...
                handles.Freq1 handles.Freq2 handles.Period handles.Calib ...
                handles.chan1 handles.chan2 handles.cal], 'enable','off')
            set([handles.mlsorder handles.nraver handles.nrreps ...
                handles.mlsordertext handles.nravertext handles.nrrepstext],...
                'visible','on'); 
            set([handles.toneperiod handles.toneperiodtext
                handles.EQfilt handles.EQfilttext],...
                'visible','off');
            set(handles.Start, 'CallBack', {@MLSmeasurement},'string','MLS GO!')
            x = 0;
    end
    handles.Data.OutputSig = x;
    guidata(Obj,handles)
end

function MLSmeasurement(Obj,Event)
    handles = guidata(Obj);
    load EQFilt.mat
    set(Obj,'enable','off')
    fs = handles.Data.FS;
    mlsorder = str2double(get(handles.mlsorder, 'string'));
    nraver = str2double(get(handles.nraver, 'string'));
    nrreps = str2double(get(handles.nrreps, 'string'));
    disp('running playrecmls for output 1 input 1')
    h(:,1) = playrecmls(fs,mlsorder,nraver,nrreps,1,1,[]);
    disp('running playrecmls for output 1 input 2')
    h(:,2) = playrecmls(fs,mlsorder,nraver,nrreps,2,1,[]);
    disp('running playrecmls for output 2 input 1')
    h(:,3) = playrecmls(fs,mlsorder,nraver,nrreps,1,2,[]);
    disp('running playrecmls for output 2 input 2')
    h(:,4) = playrecmls(fs,mlsorder,nraver,nrreps,2,2,[]);
    disp('Done')
    PlotMLS(h,handles)
    
    MSG = sprintf(' h(:,1) -> Rec1+B&K 4157\n h(:,2) -> Rec1+MIC\n h(:,3) -> Rec2+B&K 4157\n h(:,4) -> Rec2+MIC\n'); 
    save ER_10C_inCoupler_ImpulseResponse.mat fs h MSG
    %save Equalized_ER_10C_inCoupler_ImpulseResponse.mat fs h MSG
    set(handles.Start, 'value', 0, 'enable','on') 
end

function RecordAudio(Obj,Event)
    handles = guidata(Obj);
    if handles.Data.Chan1 && handles.Data.Chan2
        RecChan = [1 2];
    elseif handles.Data.Chan1 && ~handles.Data.Chan2
        RecChan = 1;
    elseif ~handles.Data.Chan1 && handles.Data.Chan2
        RecChan = 2;
    end
    if handles.Data.Out1 && handles.Data.Out2
        PlayChan = [1 2];
    elseif handles.Data.Out1 && ~handles.Data.Out2
        PlayChan = 1;
    elseif ~handles.Data.Out1 && handles.Data.Out2
        PlayChan = 2;
    else 
        PlayChan = [0 0];
    end
    ButtonPos = get(Obj,'value');
    if ButtonPos 
        AudioData = [];
        AutoPlayrecInit
        while ButtonPos
            set(handles.Period,'enable','off')
            set(handles.Calib,'enable','off')
            set(handles.chan1,'enable','off')
            set(handles.chan2,'enable','off')
            RecLength = floor(handles.Data.T*handles.Data.FS);
            Out = zeros(RecLength,2);
            if sum(PlayChan) == 0;
                out = playrec('rec',RecLength,RecChan);
            else
                out = playrec('playrec',handles.Data.OutputSig, ...
                    PlayChan, RecLength, RecChan);
            end
            playrec('block',out);
            Out(:,RecChan) = playrec('getRec',out);

            AudioData = [AudioData; Out]; 
            PlotAudio(AudioData,handles);
            ButtonPos = get(Obj,'value');
        end
        %handles.Data.AudioData = AudioData;
    else
        playrec('reset')
        set(handles.Period,'enable','on')
        set(handles.Calib,'enable','on')
        set(handles.chan1,'enable','on')
        set(handles.chan2,'enable','on')
        handles.Data.AudioData = [];
    end
    
    guidata(Obj,handles)
end


function Calibrate(Obj,Event)
    handles = guidata(Obj);
    if handles.Data.Chan1 && ~handles.Data.Chan2
        RecChan = 1;
    elseif ~handles.Data.Chan1 && handles.Data.Chan2
        RecChan = 2;
    end
    AutoPlayrecInit
    RecLength = floor(2*handles.Data.FS);
    out = playrec('rec',RecLength,RecChan);
    playrec('block',out);
    y = playrec('getRec',out);
    playrec('reset')
    
    RMS = sqrt(sum(y.^2)/length(y));
    y_dB = 10*log10(RMS^2/4e-10);
    Corr = 10^((handles.Data.CaldBSelect - y_dB(1))/20);
    y = y.*Corr; %Time signal in Pa. 
    
   
    RMS_ = sqrt(sum(y.^2)/length(y));
    y_dB_ = 10*log10(RMS_^2/4e-10);
    
    switch RecChan
        case 1
            CalInfo = sprintf('Correction Factor: %4.5f \n Measured SPL = %4.1f dB \n Adjusted SPL %4.1f', Corr, y_dB, y_dB_);
            set(handles.caltext1,'string',CalInfo);
            handles.Data.Correction(1) = Corr;
        case 2
            CalInfo = sprintf('Correction Factor: %4.5f \n Measured SPL = %4.1f dB \n Adjusted SPL %4.1f', Corr, y_dB, y_dB_);
            set(handles.caltext2,'string',CalInfo);
            handles.Data.Correction(2) = Corr;
    end         
    
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
% clear graphs:
set([handles.l1ax1 handles.l2ax1 handles.l3ax1 handles.l4ax1,...
    handles.l1ax2 handles.l2ax2 handles.l3ax2 handles.l4ax2],'visible','off')
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
legend(ax1,'off')
legend(ax2,'off')
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

y(:,1) = y(:,1).*Correction(1); %Time signal in Pa. 
y(:,2) = y(:,2).*Correction(2); %Time signal in Pa. 

% Calculate RMS and dB.
RMS_ = sqrt(sum(y.^2)/NrSamples); 
y_dB_ = 10*log10(RMS_.^2/4e-10);
y_max = max(abs(y));
y_dB_max = 10*log10(y_max.^2/4e-10);
%%
% plot in the time domain
%plot(ax1,t,y,[t(1) t(end)],[RMS_; RMS_],'linewidth',2);
set(ax1, 'xlim',[0 T],'ylim',[-max(max(abs(y))), max(max(abs(y)))])
set(handles.l1ax1,'xdata', t, 'ydata', y(:,1),'linewidth', 2,'visible','on')
set(handles.l2ax1,'xdata', t, 'ydata', y(:,2),'linewidth', 2,'visible','on')
%%
% Make the frequency analysis with the fft fuction (Fast Fourier Transform)
% and create a frequency vector for plotting.
SCDelay = floor(0.04*FS); % sound card delay 
if get(handles.stim,'value') == 1 % tone length 
    TP = floor(str2double(get(handles.toneperiod, 'string'))*FS);
    y = y(SCDelay:TP+SCDelay,:);
elseif get(handles.stim,'value') == 2 % pulse length
    TP = floor(0.02*FS);
    y = y(SCDelay:TP+SCDelay,:);
else        % only input. 
    TP = NrSamples;
end;
Y = 10*log10((2.*abs(fft(y))./length(y)).^2./4e-10);
F = 0:FS/length(y):FS-FS/length(y);
%%
% plot the frequency amplitude in dB re. 20 muPa.  
%semilogx(ax2, F(1:round(end/2)),Y(1:round(end/2),:))
set(ax2, 'ylim',[0,110])
set(handles.l1ax2,'xdata', F(1:round(end/2)), 'ydata', Y(1:round(end/2),1),'visible','on');
set(handles.l2ax2,'xdata', F(1:round(end/2)), 'ydata', Y(1:round(end/2),2),'visible','on');

%%
% Write on the figure how many seconds have been recorded:
RecSecs = TotSamples/FS;
AudioInfo = sprintf('Total Running Time: %6.2f s',RecSecs); 
set(handles.runtime,'string',AudioInfo);

IN1 = sprintf('Input 1: \n RMS: %4.2f [Pa]\n SPL: %4.1f [dB]\n SPL peak: %4.1f',...
    RMS_(1),y_dB_(1),y_dB_max(1));
IN2 = sprintf('Input 2: \n RMS: %4.2f [Pa]\n SPL: %4.1f [dB]\n SPL peak: %4.1f',...
    RMS_(2),y_dB_(2),y_dB_max(2));
set(handles.info1,'string',IN1);
set(handles.info2,'string',IN2);
drawnow
end

function PlotMLS(h,handles)
% 
FS = handles.Data.FS;
TotSamples = length(h);

ax1 = handles.Ax1;
ax2 = handles.Ax2;

t = 1/FS:1/FS:TotSamples/FS;

% convert to pascals
%Correction = handles.Data.Correction;

%y(:,1) = y(:,1).*Correction(1); %Time signal in Pa. 
%y(:,2) = y(:,2).*Correction(2); %Time signal in Pa. 

% Calculate RMS and dB.
% RMS_ = sqrt(sum(y.^2)/NrSamples); 
% y_dB_ = 10*log10(RMS_.^2/4e-10);
% y_max = max(abs(y));
% y_dB_max = 10*log10(y_max.^2/4e-10);
%%
% plot in the time domain
[Max I1d] = max(abs(h));
[Max_ I2d] = max(max(abs(h)));
premax = 20;
resplength = 2000;
set(ax1,'xlim',[0 t(end)], 'ylim',[-max(max(abs(h))), max(max(abs(h)))]);%'xlim',[t(I1d(I2d)-premax) t(I1d(I2d)+resplength-premax)],
set(handles.l1ax1,'xdata', t, 'ydata', h(:,1),'linewidth', 2,'visible','on')
set(handles.l2ax1,'xdata', t, 'ydata', h(:,2),'linewidth', 2,'visible','on')
set(handles.l3ax1,'xdata', t, 'ydata', h(:,3),'linewidth', 2,'visible','on')
set(handles.l4ax1,'xdata', t, 'ydata', h(:,4),'linewidth', 2,'visible','on')

legend(ax1,'Out1 In1','Out1 In2', 'Out2 In1', 'Out2 In2')
%%
% Make the frequency analysis with the fft fuction (Fast Fourier Transform)
% and create a frequency vector for plotting.
y = h;
nfft = length(y); %4800;
%y(:, [1 3]) = h(I1d(I2d)-premax:I1d(I2d)+resplength-premax,[1 3]);
%y(:,[2 4]) = h(I1d(2)-premax:I1d(2)+resplength-premax,[2 4]);

Y = 20*log10(abs(fft(y,nfft)));
Y(:,1) = Y(:,1)-Y(101,1);
Y(:,2) = Y(:,2)-Y(101,2);
Y(:,3) = Y(:,3)-Y(101,3);
Y(:,4) = Y(:,4)-Y(101,4);
F = 0:FS/nfft:FS-FS/nfft;
%%
% plot the frequency amplitude in dB re. 20 muPa.  
%semilogx(ax2, F(1:round(end/2)),Y(1:round(end/2),:))
set(ax2, 'ylim',[-20 20],'xgrid','on','ygrid','on')
set(handles.l1ax2,'xdata', F(1:round(end/2)), ...
    'ydata', Y(1:round(end/2),1),'visible','on');
set(handles.l2ax2,'xdata', F(1:round(end/2)), ...
    'ydata', Y(1:round(end/2),2),'visible','on');
set(handles.l3ax2,'xdata', F(1:round(end/2)), ...
    'ydata', Y(1:round(end/2),3),'visible','on');
set(handles.l4ax2,'xdata', F(1:round(end/2)), ...
    'ydata', Y(1:round(end/2),4),'visible','on');
legend(ax2,'Out1 In1','Out1 In2', 'Out2 In1', 'Out2 In2','location','southwest')
%%
% Write on the figure how many seconds have been recorded:
% RecSecs = TotSamples/FS;
% AudioInfo = sprintf('Total Running Time: %6.2f s',RecSecs); 
% set(handles.runtime,'string',AudioInfo);
% 
% IN1 = sprintf('Input 1: \n RMS: %4.2f [Pa]\n SPL: %4.1f [dB]\n SPL peak: %4.1f',...
%     RMS_(1),y_dB_(1),y_dB_max(1));
% IN2 = sprintf('Input 2: \n RMS: %4.2f [Pa]\n SPL: %4.1f [dB]\n SPL peak: %4.1f',...
%     RMS_(2),y_dB_(2),y_dB_max(2));
% set(handles.info1,'string',IN1);
% set(handles.info2,'string',IN2);
% drawnow

end
