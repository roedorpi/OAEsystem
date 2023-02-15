function checkfit

% This function runs a probe check-fit procedure, which should be included 
% before any OAE measurement. It sends an impulse and displays response of
% the system (in time and frequency domain). Output of this function
% (check-fit response and check-fit response spectrum) can be used for
% check-fit monitoring in TEOAE measurement function.
% This function also runs AutoPlayrecInit function, which initialises 
% playrec with Edirol UA-25EX sound card.

global fs playChanList recChanList check_level check_l delay_samples...
    checkFigure check_response check_impulse resp_spec chwindow scROUT...
    scLOUT scIN mic_sens speakerL speakerR PatientID Ear check_level_adj

check_w=80e-6;


%% load previous ear if any
HearingHistory = get_patientinfo(PatientID,'../emails/');
if strcmp(Ear,'L')
    if isfield(HearingHistory,'ProbeFitL')
        FittingNr = size(HearingHistory.ProbeFitL,2) + 1;
    else
        FittingNr = 1;
    end
elseif strcmp(Ear,'R')
    if isfield(HearingHistory,'ProbeFitR')
        FittingNr = size(HearingHistory.ProbeFitR,2) + 1;
    else
        FittingNr = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if playChanList == 1; scOUT = scLOUT; else scOUT = scROUT; end
if playChanList == 1; speak_sens = speakerL; else speak_sens = speakerR; end

AutoPlayrecInit(48e3);

%% check fit

check_Pa=20*10^((check_level/20)-6);
check_digi=check_Pa/(speak_sens*scOUT);

check_impulse=ones(1,round(check_w*fs));
check_impulse=[zeros(1,100) check_impulse zeros(1,round((check_l-check_w)*fs)-100)]; %impulse length 0.1 s

%compensating for the loudspeaker transfer function
check_impulse = lsinvopt10kHz(check_impulse')';

% low pass filtering
load('lowfilter20db.mat')
check_impulse=filter(hd,check_impulse);

% adjust the level after filtering
check_impulse_peak = max(abs(check_impulse));
correction = check_digi/check_impulse_peak;
check_impulse=check_impulse.*correction;

% initialize checkfit figure
scrsize = get(0,'ScreenSize');
checkFigure=figure;
set(checkFigure,'OuterPosition',[0 30 scrsize(3) scrsize(4)-30])
setappdata(checkFigure,'stopp',0);

t=(0:1/fs:(length(check_impulse)-1)/fs)*1000;

TimePlot = axes('parent',checkFigure,...
    'OuterPosition',[0 0.5 0.7 0.5],...
    'XLim',[0 100], 'Ylim', [-0.2, 0.2],'nextplot','add');
ylabel('Checkfit response [Pa]')
xlabel('Time [ms]')

checkFFTAxes = axes('parent',checkFigure,...
    'OuterPosition',[0 0 0.7 0.5],'nextplot','add');

checktekst1 = uicontrol(checkFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.6 0.25 0.12],...
    'BackgroundColor','w',...
    'Fontsize',11);

startbut = uicontrol(checkFigure,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.7 0.4 0.25 0.1],...
    'BackgroundColor','w',...
    'Fontsize',11,...
    'String','Press to start the measurement',...
    'Callb',@(obj,evnt)stopper(obj,evnt));

uicontrol(startbut) % give start button focus

if FittingNr > 1
    if strcmp(Ear,'L')
        plot(TimePlot,t,HearingHistory.ProbeFitL(:,end),...
            'color',[0.8 0.8 0.8],'linewidth',2);
        org_resp_peak=20*log10(max(abs(HearingHistory.ProbeFitL(:,end)))/20e-6);
        nfft_=length(HearingHistory.ProbeFitL(:,end));
        f_=0:fs/nfft_:(fs/2-fs/nfft_);
        org_spec=abs(fft(HearingHistory.ProbeFitL(:,end),nfft_))/(1e-3*fs);
        org_spec=2*org_spec(1:nfft_/2);
        plot(checkFFTAxes,f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
            'linewidth',3,'color',[0.7 0.7 0.7]);
    elseif strcmp(Ear,'R')
        plot(TimePlot,t,HearingHistory.ProbeFitR(:,end),...
            'color',[0.8 0.8 0.8],'linewidth',2);
        org_resp_peak=20*log10(max(abs(HearingHistory.ProbeFitR(:,end)))/20e-6);
        nfft_=length(HearingHistory.ProbeFitR(:,end));
        f_=0:fs/nfft_:(fs/2-fs/nfft_);
        org_spec=abs(fft(HearingHistory.ProbeFitR(:,end),nfft_))/(1e-3*fs);
        org_spec=2*org_spec(1:nfft_/2);
        plot(checkFFTAxes,f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
            'linewidth',3,'color',[0.7 0.7 0.7]);
    end
else
    org_resp_peak = 0;
end
checkWaveLine = line('XData',t,'YData',zeros(1,length(t)),'parent',TimePlot);
hSpecArea = area(checkFFTAxes,[1 2], [3 4]);
check_level_adj = check_level;
while getappdata(checkFigure,'stopp')==0
    pageNumber = playrec('playrec',check_impulse',playChanList,length(check_impulse),recChanList);
    playrec('block',pageNumber);
    check_response = playrec('getRec',pageNumber);
%    check_response = hpf(check_response,500,fs);
    win=hann(400);
    chwindow=[win(1:length(win)/2); ones(length(check_response)-length(win),1);...
    win(length(win)/2+1:length(win))];
    check_response=check_response.*chwindow;    
    
    
    
    check_response=check_response/(scIN*mic_sens); %[Pa]
    check_resp_peak=20*log10(max(abs(check_response))/20e-6);
    
  
    
    nfft=length(check_response);
    f=0:fs/nfft:(fs/2-fs/nfft);
    resp_spec=abs(fft(check_response',nfft))/(1e-3*fs);
    resp_spec=2*resp_spec(1:nfft/2);

    chpeak=sprintf(['\nCheckfit response peak level =\n'...
        num2str(round(10*check_resp_peak)/10)...
            ' dB SPL \n'...
            'First fitting response peak level = \n'...
         num2str(round(10*org_resp_peak)/10)...
            ' dB SPL']);
    set(checktekst1,'String',chpeak)
    set(checkWaveLine,'YData',check_response);
    delete(hSpecArea)
    hSpecArea = area(checkFFTAxes,f,[resp_spec(1) 20*log10(resp_spec(2:end)/20e-6)]);
    set(checkFFTAxes,...
        'xlim', [0 6000],...
        'XTick',[1000 2000 3000 4000 5000],...
        'ylim', [0 100],...
        'YTick',0:20:100,...
        'XGrid','on');
    xlabel('Frequency [Hz]')
    ylabel('Checkfit response spectrum [dBSPL]')
    grid on
    set(gca,'Layer','top')
    drawnow
    if abs(check_level - check_resp_peak) > 0.2
        dBdiff = check_level - check_resp_peak;
        
        check_Pa_new=20*10^(((check_level_adj + dBdiff)/20)-6);
        check_digi_new=check_Pa_new/(speak_sens*scOUT);
        
        check_impulse_peak_ = max(abs(check_impulse));
        correction_ = check_digi_new/check_impulse_peak_;
        check_impulse=check_impulse.*correction_;
        check_level_adj = check_level_adj + dBdiff;
    end
end
close

if strcmp(Ear,'L')
   HearingHistory.ProbeFitL(:,FittingNr) = check_response;
elseif strcmp(Ear,'R')
   HearingHistory.ProbeFitR(:,FittingNr) = check_response;
end
Name2Save = strcat('../emails/', HearingHistory.ID, '.mat');
save(Name2Save, 'HearingHistory')

% Calculate delay from checkfit data
[correlation lags] = xcorr(check_response,check_impulse);
[CorMax Idx] = max(correlation); %#ok
delay_samples = lags(Idx);

% 
% [peak peak_index] = max(check_response); %#ok
% 
% k1=1;
% n=20;
% M=zeros(1,floor(length(check_response)/n));
% for i=1:n:length(check_response)-n
%     if k1==floor(length(check_response)/n)
%         M(k1)=mean(abs(check_response(i:end)));
%     else  M(k1)=mean(abs(check_response(i:i+n)));
%     end
%     k1=k1+1;
% end
% p=round(peak_index/n);
% k2=p;
% while abs(M(k2)-M(k2-1)) > 0.05*M(p)
%     k2=k2-1;
% end
% delay_samples=k2*n-100;
% if delay_samples<0, delay_samples=0; end

end



function stopper(varargin)
    setappdata(get(varargin{1},'parent'),'stopp',1);
    
end
