function delay_samples = appcheckfit(app)

% This function runs a probe check-fit procedure, which should be included 
% before any OAE measurement. It sends an impulse and displays response of
% the system (in time and frequency domain). Output of this function
% (check-fit response and check-fit response spectrum) can be used for
% check-fit monitoring in TEOAE measurement function.



%% load previous ear if any
if strcmp(app.Data.Ear,'L')
    if isfield(app.Data,'ProbeFitL')
        FittingNr = size(app.Data.ProbeFitL,2) + 1;
    else
        FittingNr = 1;
    end
elseif strcmp(Ear,'R')
    if isfield(app.Data,'ProbeFitR')
        FittingNr = size(app.Data.ProbeFitR,2) + 1;
    else
        FittingNr = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scDAC = 1; 
fs = app.AuIO.SampleRate;
bf = app.AuIO.BufferSize;    
check_w = 80e-6;
check_l = app.Config.CheckFit.Impulselength;
check_level = app.Config.CheckFit.dBSPL;
  
%% check fit

check_Pa=20*10^((check_level/20)-6);
check_digi=check_Pa/scDAC;

check_impulse=ones(1,round(check_w*fs));
check_impulse=[zeros(1,100) check_impulse zeros(1,round((check_l-check_w)*fs)-100)]; 
%check_impulse = check_impulse(1:bf);
% %compensating for the loudspeaker transfer function
% check_impulse = lsinvopt10kHz(check_impulse')';

% low pass filtering
load('lowfilter20db.mat')
check_impulse=filter(hd,check_impulse);

% adjust the level after filtering
check_impulse_peak = max(abs(check_impulse));
correction = check_digi/check_impulse_peak;
check_impulse=check_impulse.*correction;
x = playRecSig(app.AuIO,repmat(check_impulse',1,length(app.AuIO.PlayerChannelMapping)));
t=(0:1/fs:(length(x)-1)/fs)*1000;

if FittingNr > 1
    if strcmp(app.Data.Ear,'L')
        t_=(0:1/fs:(length(app.Data.ProbeFitL)-1)/fs)*1000; 
        plot(app.UIAxesCheckFitTime,t_,app.Data.ProbeFitL(:,end),...
            'color',[0.8 0.8 0.8],'linewidth',2);
        org_resp_peak=20*log10(max(abs(app.Data.ProbeFitL(:,end)))/20e-6);
        nfft_=length(app.Data.ProbeFitL(:,end));
        f_=0:fs/nfft_:(fs/2-fs/nfft_);
        org_spec=abs(fft(app.Data.ProbeFitL(:,end),nfft_))/(1e-3*fs);
        org_spec=2*org_spec(1:nfft_/2);
        plot(app.UIAxesCheckFitFreq,f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
            'linewidth',3,'color',[0.7 0.7 0.7]);
    elseif strcmp(app.Data.Ear,'R')
        plot(app.UIAxesCheckFitTime,t,app.Data.ProbeFitR(:,end),...
            'color',[0.8 0.8 0.8],'linewidth',2);
        org_resp_peak=20*log10(max(abs(app.Data.ProbeFitR(:,end)))/20e-6);
        nfft_=length(app.Data.ProbeFitL(:,end));
        f_=0:fs/nfft_:(fs/2-fs/nfft_);
        org_spec=abs(fft(app.Data.ProbeFitR(:,end),nfft_))/(1e-3*fs);
        org_spec=2*org_spec(1:nfft_/2);
        plot(app.UIAxesCheckFitFreq,f_,[org_spec(1) 20*log10(org_spec(2:end)'/20e-6)],...
            'linewidth',3,'color',[0.7 0.7 0.7]);
    end
else
    org_resp_peak = 0;
end
checkWaveLine = line('XData',t,'YData',zeros(1,length(t)),'parent',app.UIAxesCheckFitTime);
hSpecArea = area(app.UIAxesCheckFitFreq,[1 2], [3 4]);
check_level_adj = check_level;
mouse_pressed = false;
while ~mouse_pressed
    MoPos = app.UIAxesCheckFitTime.CurrentPoint;
    if MoPos(1,1) >  app.UIAxesCheckFitTime.XLim(1) &&  ...
            MoPos(1,1) <  app.UIAxesCheckFitTime.XLim(2) && ...
            MoPos(1,2) >  app.UIAxesCheckFitTime.YLim(1) &&  ...
            MoPos(1,2) <  app.UIAxesCheckFitTime.YLim(2)
        mouse_pressed = true;
    end

    check_response = playRecSig(app.AuIO,[check_impulse' zeros(length(check_impulse),1)]);
    check_response = hpf(check_response(:,1),500,fs);
    win=hann(400);
    chwindow=[win(1:length(win)/2); ones(length(check_response)-length(win),1);...
    win(length(win)/2+1:length(win))];

%     [correlation lags] = xcorr(check_response,check_impulse);
%     [CorMax Idx] = max(correlation); %#ok
%     delay_samples = lags(Idx)/fs;

    check_response=check_response.*chwindow;    
    check_response=check_response/scDAC; %[Pa]
    check_resp_peak=20*log10(max(abs(check_response))/20e-6);
    
    nfft=length(check_response);
    f=0:fs/nfft:(fs/2-fs/nfft);
    resp_spec=abs(fft(check_response',nfft))/(1e-3*fs);
    resp_spec=2*resp_spec(1:nfft/2);
    
    app.MeasInfoTab.Data(6,1) = {roundn(check_resp_peak,-1)};
    app.MeasInfoTab.Data(7,1) = {roundn(10*org_resp_peak,1)};
    set(checkWaveLine,'YData',check_response);
    set(app.UIAxesCheckFitTime,...
        'xlim', [0 1.5*check_l*1000],...
        'ylim', [-0.2 0.2],...
        'XGrid','on');
    delete(hSpecArea)
    hSpecArea = area(app.UIAxesCheckFitFreq,f,[resp_spec(1) 20*log10(resp_spec(2:end)/20e-6)]);
    set(app.UIAxesCheckFitFreq,...
        'xlim', [0 6000],...
        'XTick',[1000 2000 3000 4000 5000],...
        'ylim', [0 100],...
        'YTick',0:20:100,...
        'XGrid','on');
    set(gca,'Layer','top')
    drawnow
    if abs(check_level - check_resp_peak) > 0.2
        dBdiff = check_level - check_resp_peak;
        
        check_Pa_new=20*10^(((check_level_adj + dBdiff)/20)-6);
        check_digi_new=check_Pa_new/scDAC;
        
        check_impulse_peak_ = max(abs(check_impulse));
        correction_ = check_digi_new/check_impulse_peak_;
        check_impulse=check_impulse.*correction_;
        check_level_adj = check_level_adj + dBdiff;
    end
end

if strcmp(app.Data.Ear,'L')
   app.Data.ProbeFitL(:,FittingNr) = check_response;
elseif strcmp(app.Data.Ear,'R')
   app.Data.ProbeFitR(:,FittingNr) = check_response;
end

% Calculate delay from checkfit data
[correlation lags] = xcorr(check_response,check_impulse);
[CorMax Idx] = max(correlation); %#ok
delay_samples = lags(Idx)/fs;
