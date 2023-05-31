function TEOAEmeas(varargin)
%% Main figure:

scrsize = get(0,'ScreenSize');

main = figure(...
    "OuterPosition",[0 60 scrsize(3)/2 scrsize(4)-60],...
    'ToolBar','none','MenuBar','none');%,'DeleteFcn', @CleanUp);
for i = 1:4
    ax(i) = axes('Parent',main,'Units','normalized','box','on');
end

%%
ax(1).Position = [0.25 0.54 0.49 0.38];
ax(2).Position = [0.25 0.08 0.49 0.38];
ax(3).Position = [0.775 0.29 0.2 0.18];
ax(4).Position = [0.775 0.07 0.2 0.18];
ax(3).ButtonDownFcn = @StopCheckFit;
ax(1).Title.String = 'Time Signal';
ax(1).XLabel.String = 'Time [ms]';
ax(1).YLabel.String = 'Sound Pressure [Pa]';
ax(2).Title.String = 'Frequency Response';
ax(2).XLabel.String = 'Frequency [kHz]';
ax(2).YLabel.String = 'Amplitude [dB]';
ax(3).Title.String = 'Probe Fit';
ax(3).XLabel.String = 'Time [s]';
ax(3).YLabel.String = 'Sound Pressure [Pa]';
ax(4).XLabel.String = 'Frequency [kHz]';
ax(4).YLabel.String = 'Amplitude [dB]';

set(ax(2),'ylim',[-60 10], 'xlim', [0 6000],'xtick',1000:1000:6000,...
    'xticklabel',{'1' '2' '3' '4' '5' ''},'xgrid','on','ygrid','on',...
    'box','on');
set(ax(1),'ylim',[-0.05 0.05], 'xlim', [0 20],'xgrid','on','ygrid','on',...
    'box','on');

%% measurement information
MeasInfoTab = uitable(main,'Units','normalized'); 
MeasInfoTab.Position = [0.775 0.52 0.2 0.4];
MeasInfoTab.ColumnName = {'Parameter','Value'};
MeasInfoTab.ColumnWidth = {180, 120};

tabData = {'Date:'; 'Time:'; 'TEOAE Level [dB]:';...
    'Noise Level [dB]:'; 'Percent Rejected [%]:'; ...
    'Wave Repro [%]:'; 'AB Correlation [%]';...
    'Check Fit Peak Level [dB]:';...
    'Prev. Check Fit Peak Level [dB]:';'Check Fit Correlation:';...
    'Monitor/Original check fit peak [dB]'};

MeasInfoTab.Data = cell(length(tabData),2);
MeasInfoTab.Data(:,1) = tabData; 

MeasInfoTab.Data(1,2) = cellstr(datestr(floor(now)));


%% Measurement parameters
MeasParamTab = uitable(main,'Units','normalized','ColumnEditable',logical([0,1]));
MeasParamTab.Position = [0.02 0.35 0.19 0.57];
MeasParamTab.ColumnName = {'Parameter','Value'};
MeasParamTab.ColumnWidth = {180, 100};
MeasParamTab.Data = cell(11,2);
MeasParamTab.Data(:,1) = {'Max peak Level [dB SPL]';...
    'Noise rejection level [dB SPL]'; 'Time rejection [ms]';...
    'Stimulus spacing [ms]'; 'Number of measurements'; ...
    'High-pass cutoff [Hz]'; 'Low-pass cutoff [Hz]';...
    'Checkfit impulse level [dB SPL]';'Checkfit monitor frequency [Hz]';...
    'Contralateral Stimulation'; 'Ear'};
% Default Values
MeasParamTab.Data(:,2) = {82;40;4;20;512;500;10000;80;20;false;'L'};
%MeasParamTab.CellEditCallback %= @updateMeasParam;     
%%
ButtonControl = uicontrol(main,"Style","pushbutton","Units","normalized",...
    "String","Start","Callback",@runcheckfit);
ButtonControl.Position = [0.02 0.07 0.19 0.1];

%%
NewSubjectButton = uicontrol(main,"Style","pushbutton","Units","normalized",...
    "String","New Subject","Callback",@GetSetSubjectFile,'tag','NewSubject');
NewSubjectButton.Position = [0.02 0.2 0.09 0.1];
LoadSubjectButton = uicontrol(main,"Style","pushbutton","Units","normalized",...
    "String","Load Subject","Callback",@GetSetSubjectFile,'tag','LoadSubject');
LoadSubjectButton.Position = [0.12 0.2 0.09 0.1];

%% Data structure
Data = struct(...
    'A',[],'B',[],'Stimulus',[],'Ear','','PatientID','',...
    'Date','','StartTime','','EndTime','','ProbeFitL',[],...
    'ProbeFitR',[]);
Config = struct(...
    'Fs', [],'TEdBPeak', [],'NumofMeasurements', [],...
    'TimeRejection', [],'NoiseThreshold', [],...
    'BandPassFilter', [], 'CheckFit', struct,'TEsignal',[],'CLS',[]);
Config.CheckFit = struct('scDAC', 1,'scADC', 0.5,'scDelay',[],...
    'check_w', 80e-6,'check_l', 0.1);

%% Calibration 
load('ER10Bch1_KemarrightEarch2.mat','ChCalRMS');
load("CalibFile_1-5-23_ER10Bch1_KemarchRightch2_IP30right.mat","CalibData");
% find the value closes to 94 dBSPL at 1000kHz in the Right channel of the
% CalibData table. When the earphone is actuated with a signal of -XX dBFS
% we get XX dB SPL measured at the couplers microphone (Right Channel in 
% table). 
FSLev = -18;
SPLLev = 93.1356;
Pa = 10^(SPLLev/20)*20e-6;
Amp = 10^(FSLev/20);
Config.CheckFit.scDAC = Pa/Amp;
Config.CheckFit.scADC = ChCalRMS(1); %Channel 1 is the probe mic. 

    
%% Generate figure data

h.Data = Data;
h.Config = Config;
h.main = main;
h.ax = ax;
h.ButtonControl = ButtonControl;
h.MeasInfoTab = MeasInfoTab;
h.MeasParamTab = MeasParamTab;

setappdata(main,'main',h)
AudioPlayerInit(main)
end

function GetSetSubjectFile(obj,event)
%% Open or make new response file
h = getappdata(obj.Parent,'main');

if strcmp(obj.Tag,'LoadSubject')
    [File.Name, File.Path] = uigetfile('data\*.mat','Select existing data file.');
else
    [File.Name, File.Path] = uiputfile('data\*.mat','Type the name of a new data file.');
end
if File.Name ~= 0
    FileIO = matfile([File.Path,File.Name],"Writable",true);
 else
    uiwait(msgbox(main,'No Patient file was loaded or created, try again.'))
    FileIO = 0;

end
h.File = File;
h.FileIO = FileIO;

setappdata(obj.Parent,'main',h);
end

function AudioPlayerInit(obj)
%% Initialization of full duplex sound card using ASIO driver
% duplex audio object
h = getappdata(obj,'main');
AuIO = audioPlayerRecorder;
b = AuIO.getAudioDevices;
AuIO.Device = b{2};
AuIO.SampleRate = 48000;
AuIO.PlayerChannelMapping = [1 2];
AuIO.RecorderChannelMapping = [1 2];
% Check the buffer size properties on the ASIO settings to set the same
% here.
AuIO.SupportVariableSize = true;
AuIO.BufferSize = 1024;
% play zeros to lock device for use. 
try
    [~] = AuIO(0.6*zeros(10*AuIO.BufferSize,length(AuIO.PlayerChannelMapping)));
catch
    sprintf('Device: %s is locked or cannot be initialized',AuIO.Device)
    
    %uiconfirm(h.main,MsG,'','Options',{'Ok.'},'icon','error')
end
h.AuIO = AuIO;
setappdata(obj,'main',h);

end

function runcheckfit(obj,event)
    fig = gcbf;
    h = getappdata(fig,'main');
    if strcmp(h.ButtonControl.String,'Start')
        h.ButtonControl.String = 'Stop';
        h.ButtonControl.BackgroundColor = [0.8 0.2 0.2];
    
        h.MeasInfoTab.Data(2,2) = cellstr(datestr(rem(now,1),'HH:MM:SS'));
        fs = h.AuIO.SampleRate;
        delete(h.ax(3).Children);
        h.ax(3).NextPlot = 'add';
        delete(h.ax(4).Children);
        h.ax(4).NextPlot = 'add';
        %% load previous ear if any
        if strcmp(h.MeasParamTab.Data{11,2},'L')
            FittingNr = size(h.Data.ProbeFitL,2) + 1;
        else
            FittingNr = size(h.Data.ProbeFitR,2) + 1;
        end
        if FittingNr > 1
            if strcmp(h.MeasParamTab.Data{11,2},'L')
                t_=(0:1/fs:(length(h.Data.ProbeFitL)-1)/fs)*1000; 
                nfft_=length(h.Data.ProbeFitL(:,end));
                f_=0:fs/nfft_:(fs/2-fs/nfft_);
                plot(h.ax(3),t_,h.Data.ProbeFitL(:,end),...
                    'color',[0.8 0.8 0.8],'linewidth',2);
                org_resp_peak = 20*log10(max(abs(h.Data.ProbeFitL(:,end)))/20e-6);
                
                org_spec=abs(fft(h.Data.ProbeFitL(:,end),nfft_))/length(h.Data.ProbeFitL);
                org_spec=2*org_spec(1:nfft_/2);
                plot(h.ax(4),f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
                    'linewidth',3,'color',[0.7 0.7 0.7]);
            elseif strcmp(h.MeasParamTab.Data{11,2},'R')
                t_=(0:1/fs:(length(h.Data.ProbeFitR)-1)/fs)*1000; 
                nfft_=length(h.Data.ProbeFitR(:,end));
                f_=0:fs/nfft_:(fs/2-fs/nfft_);
        
                plot(h.ax(3),t_,h.Data.ProbeFitR(:,end),...
                    'color',[0.8 0.8 0.8],'linewidth',2);
                org_resp_peak=20*log10(max(abs(h.Data.ProbeFitR(:,end)))/20e-6);
                
                org_spec=abs(fft(h.Data.ProbeFitR(:,end),nfft_))/length(h.Data.ProbeFitR);
                org_spec=2*org_spec(1:nfft_/2);
                plot(h.ax(4),f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
                    'linewidth',3,'color',[0.7 0.7 0.7]);
            end
        else
            org_resp_peak = 0;
        end
        
        h.MeasInfoTab.Data(9,2) = {roundn(org_resp_peak,1)};
        %%
        % Output and Input sencitivity needs calibration
        scDAC = h.Config.CheckFit.scDAC;
        scADC = h.Config.CheckFit.scADC;
        %stimulus
        check_w = h.Config.CheckFit.check_w;
        check_level = h.MeasParamTab.Data{8,2};
        check_l = h.Config.CheckFit.check_l; 
        
        check_Pa=20*10^((check_level/20)-6);
        check_digi=check_Pa/scDAC;
        
        check_impulse=ones(1,round(check_w*fs));
        check_impulse=[zeros(1,100) check_impulse zeros(1,round((check_l-check_w)*fs)-100)]; %impulse length 0.1 s
        
        % compensating for the loudspeaker transfer function
        check_impulse = lsinvopt10kHz(check_impulse')';
        
        % low pass filtering
        load('lowfilter20db.mat')
        check_impulse=filter(hd,check_impulse);
        
        % adjust the level after filtering
        check_impulse_peak = max(abs(check_impulse));
        correction = check_digi/check_impulse_peak;
        check_impulse=check_impulse.*correction;
        
        h.Config.CheckFit.check_impulse = check_impulse;
        %% 
        x = playRecSig(h.AuIO,zeros(length(check_impulse),2),0);
        t = (1/fs:1/fs:length(x)/fs)*1000;
        %plot only channel 1
        pl = line('XData',t,'YData',x(:,1),'parent',h.ax(3));
        nfft_=length(t);
        f_=0:fs/nfft_:(fs/2-fs/nfft_);
        pf = area(h.ax(4),f_,40*ones(length(f_),1));
        
        %%
        h.mouse_pressed = false;
        setappdata(fig,'main',h)
        %h = getappdata(fig,'main');
        while ~h.mouse_pressed
            
            x_ = playRecSig(h.AuIO,[check_impulse',zeros(length(check_impulse),1)],0);
            scDelay = length(x_) - length(check_impulse);
            
        
            %x = mean(cat(3,x,x_),3);
            check_response = hpf(x_(:,1),500,fs);
            win=hann(400);
            chwindow=[win(1:length(win)/2); ones(length(check_response)-length(win),1);...
            win(length(win)/2+1:length(win))];
        
            check_response=check_response.*chwindow;    
            check_response=check_response/scADC; %[Pa]
            check_resp_peak=20*log10(max(abs(check_response))/20e-6);
            
            nfft=length(check_response);
            f=0:fs/nfft:(fs/2-fs/nfft);
            resp_spec=abs(fft(check_response',nfft))/(length(check_response));
            resp_spec=2*resp_spec(1:nfft/2);
            
            h.MeasInfoTab.Data(8,2) = {roundn(check_resp_peak,-1)};
            
            pl.YData = check_response;
            set(h.ax(3),...
                'xlim', [0 1.1*check_l*1000],...
                'ylim', [-0.25 0.25],...
                'XGrid','on');
            
            pf.YData = [resp_spec(1) 20*log10(resp_spec(2:end)/20e-6)];
            set(h.ax(4),...
                'xlim', [0 6000],...
                'XTick',[1000 2000 3000 4000 5000],...
                'ylim', [0 80],...
                'YTick',0:20:100,...
                'XGrid','on');
            set(gca,'Layer','top')
            setappdata(fig,'main',h);
            drawnow
            h = getappdata(fig,'main');
        end
        if strcmp(h.MeasParamTab.Data{11,2},'L') 
            h.Data.ProbeFitL = [h.Data.ProbeFitL check_response];
        else
            h.Data.ProbeFitR = [h.Data.ProbeFitR check_response];
        end
        setappdata(fig,'main',h);
        if strcmp(h.ButtonControl.String,'Stop')
            teoae_recording_new(h.main);
        end
    else
        
        h.ButtonControl.String = 'Start';
        h.ButtonControl.BackgroundColor = [0.2 0.8 0.2];
        
    end
end

function StopCheckFit(obj,~)
    f = gcbf; 
    h = getappdata(f,'main');    
    h.mouse_pressed = true;
    h.AuIO.release;
    setappdata(f,'main',h);

end