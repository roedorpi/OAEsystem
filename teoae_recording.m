function teoae_recording(~,~)

% This function masures TEOAE using non-linear mode of stimuli. Has to be 
% run after checkfit.m, because of the check-fit monitoring included in it.
% Input parameters are checkfit response, delay (from checkfit function),
% input and output channels, stimuli level, number of measurements,
% rejection time, rejection threshold, band pass filtering cutoff
% frequencies.
% The function outputs two matrices, with time signals from all the
% measurements. The TEOAE processing funcion can be run afterwards to get
% more TEOAE parameters (mean response, reproducibility, etc.).

global A B check_response resp_spec check_impulse check_level check_level_adj...
    delay_samples checkfit_monit fs playChanList recChanList ydb ...
    measurements time_rejection threshold spacing f1 f2 niceFigure monit ...
    tt ABAx AaLine BbLine scLOUT scROUT scIN speakerL speakerR mic_sens ...
    chwindow stimulus

global PatientID Ear ID

load('lowfilter20db.mat')

%%%%%%%%%%%%%%%%%%%%%%% scaling factors %%%%%%%%%%%%%%%%%%%%%%%%
% scLOUT      = 3.3220;   % [Vrms/output digital rms] output 1 (left speaker)
% scROUT      = 3.4266;   % [Vrms/output digital rms] output 2 (right speaker)
% speak_sens  = 0.07963;  % [Pa/Vrms]
% mic_sens    = 0.5;      % [Vrms/Pa]
% scIN        = 0.433;    % [digital rms/Vrms]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if playChanList == 1; scOUT = scLOUT; else scOUT = scROUT; end
if playChanList == 1; speak_sens = speakerL; else speak_sens = speakerR; end

%% Generate impulse stimuli

LevDiff = ydb - check_level;
StimMaxPeakLev =  check_level_adj + LevDiff;

yPa=20*10^((StimMaxPeakLev/20)-6);
imp_level=yPa/(speak_sens*scOUT); %[digital]

    impulse1 = zeros(1,spacing*fs);
    impulse1(1:round(80e-6*fs)) = 1/3;

    impulse3 = zeros(1,spacing*fs);
    impulse3(1:round(80e-6*fs)) = -1;

impulses = [impulse1 impulse1 impulse1 impulse3];

%compensating for the loudspeaker transfer function
impulses = lsinvopt10kHz(impulses')';

% low pass filtering
impulses=filter(hd,impulses);

% adjust the level after filtering
impulses=imp_level*impulses/max(abs(impulses)); %0.8515;

% stimulus=impulses;

%% Get the response

t=(0:1/fs:(length(check_response)-1)/fs)*1000;
nfft=length(check_response);
f=0:fs/nfft:(fs/2-fs/nfft);

%initialize the figures
scrsize = get(0,'ScreenSize');
niceFigure=figure;
set(niceFigure,'OuterPosition',[0 30 scrsize(3) scrsize(4)-30])
y_wave_limit = max(abs(check_response))+0.25*max(abs(check_response));
axes('parent',niceFigure,...
    'OuterPosition',[0.6 0.25 0.4 0.25],...
    'ylim',[-y_wave_limit y_wave_limit],'box','on');
ylabel('Amplitude [Pa]')
xlabel('Time [ms]')
line('XData',t,'YData',check_response,'Color','y');
waveLine = line('XData',t,'YData',zeros(1,length(t)));
legend('Original response','Monitoring response',...
    'Location','NorthEast')

axes('parent',niceFigure,...
    'OuterPosition',[0.6 0 0.4 0.25],...
    'xlimmode','manual',...
    'xlim',[0 6000],...
    'XTick',[1000 2000 3000 4000 5000],...
    'xticklabel',{'1' '2' '3' '4' '5' ''},...
    'ylim',[0 100],...
    'YTick',0:20:100,'box','on');
xlabel('Frequency [kHz]')
ylabel('Amplitude [dB SPL]')
line('XData',f,'YData',20*log10(resp_spec/20e-6),'Color','y');
fftLine = line('XData',f,'YData',zeros(1,length(f)));
grid on


tt=(0:1/fs:((length(impulses)/4)-1)/fs)*1000;

ABAx = axes('parent',niceFigure,...
    'OuterPosition',[0 0.5 0.65 0.5],...
    'Ylim',[-0.05 0.05],'fontsize',14,'box','on');
xlabel('Time [ms]','fontsize',14);
ylabel('Acoustic pressure [Pa]','fontsize',14);
AaLine = line('parent',ABAx,'XData',tt,'YData',zeros(1,length(tt)),'Color','b');
BbLine = line('parent',ABAx,'XData',tt,'YData',zeros(1,length(tt)),'Color','g');
legend('Averaged waveform A','Averaged waveform B','Location','NorthWest')


SPECAx = axes('parent',niceFigure,...
    'Outerposition',[0 0 0.65 0.5],'fontsize',14);

ff = (0:((length(impulses)/4))/2-1)/((length(impulses)/4)).*fs;
TESpec = area(SPECAx,ff,zeros(1,length(ff)),-40,'facecolor','g',...
    'visible','off');
ylabel(SPECAx,'SPL [dB re. 20\mu Pa]','fontsize',14);
xlabel(SPECAx,'Frequency [kHz]','fontsize',14);
set(SPECAx,'nextplot','add')
NoiSpec = area(SPECAx,ff,zeros(1,length(ff)),-40,'facecolor','r',...
    'visible','off');
set(SPECAx,'ylim',[-40 20], 'xlim', [0 6000],'xtick',1000:1000:6000,...
    'xticklabel',{'1' '2' '3' '4' '5' ''},'xgrid','on','ygrid','on',...
    'box','on');
legend(SPECAx,{'TEOAE' 'Noise'},'location','northeast','fontsize',14)
legend(SPECAx,'boxoff')


SubjInfo = ['Subj. ID: ', ID ' - Ear: ' Ear] ;
uicontrol(niceFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.9 0.25 0.05],...
    'BackgroundColor','w',...
    'Fontsize',11,...
    'String',SubjInfo,...
    'fontsize',14,...
    'fontweight','bold');

rejtekst = uicontrol(niceFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.83 0.25 0.05],...
    'BackgroundColor','w',...
    'Fontsize',11,...
    'String','Rejected responses: ');

corrtekst = uicontrol(niceFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.77 0.25 0.05],...
    'BackgroundColor','w',...
    'Fontsize',11,...
    'String','AB waveforms correlation: ');

LevText = sprintf('TEOAE level: \nNoise level:');

leveltekst = uicontrol(niceFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.71 0.25 0.05],...
    'BackgroundColor','w',...
    'Fontsize',11,...
    'String',LevText);

tekst = uicontrol(niceFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.7 0.55 0.25 0.11],...
    'BackgroundColor','w',...
    'Fontsize',11);

inboxdisp=sprintf(['Monitoring/original checkfit peak:\n' '\n'...
            'Check-fit signal correlation: ']);
set(tekst,'String',inboxdisp)

% uicontrol(niceFigure,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.72 0.6 0.22 0.05],...
%     'BackgroundColor','w',...
%     'Fontsize',11,...
%     'FontWeight','bold',...
%     'String','SAVE RESULTS',...
%     'Callback',@teoae_save);

% uicontrol(niceFigure,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.72 0.54 0.22 0.05],...
%     'BackgroundColor','w',...
%     'Fontsize',11,...
%     'FontWeight','normal',...
%     'String','BACK TO TEOAE MEASUREMENTS',...
%     'Callback',@backtoteoae);

% uicontrol(niceFigure,...
%     'Style','pushbutton',...
%     'Units','normalized',...
%     'Position',[0.72 0.48 0.22 0.05],...
%     'BackgroundColor','w',...
%     'Fontsize',11,...
%     'FontWeight','normal',...
%     'String','BACK TO MAIN MENU',...
%     'Callback',@maingui);

drawnow

k=0;
reject_count=0;
AB_corr=1;
stimulus=NaN(0.08*fs,measurements);

%filtering the response (coefficients)
f1=2*f1/fs; f2=2*f2/fs;
[coefb coefa]=butter(4,[f1 f2]);

A_sum=[]; B_sum=[];
if round(checkfit_monit/10) ~= checkfit_monit/10
    checkfit_monit = ceil(checkfit_monit/10)*10;
end

i = 0;
while k<measurements;
   %% Checkfit monitoring
   i = i + 1;
   if round((i-1)/(checkfit_monit/10))==(i-1)/(checkfit_monit/10) && i~=1
        pageNumber = playrec('playrec',check_impulse',playChanList,-1,recChanList);
        playrec('block',pageNumber);
        [monit] = playrec('getRec',pageNumber);
        
        
        %monit=hpf(monit,500,fs);
        
        win=hann(400);
        chwindow=[win(1:length(win)/2); ones(length(check_response)-length(win),1);...
        win(length(win)/2+1:length(win))];
        monit=monit.*chwindow;

        monit=monit/(mic_sens*scIN); %[Pa]
        monit_peak=20*log10(max(abs(monit))/20e-6);

        mon_spec=abs(fft(monit',nfft))/(1e-3*fs);
        mon_spec=2*mon_spec(1:length(mon_spec)/2);
        mon_check_corr=corr(monit,check_response);
        mck=['Check-fit signal correlation: ' num2str(round(100*mon_check_corr)/100)];
        check_resp_peak=20*log10(max(abs(check_response))/20e-6);
        mchpeak=sprintf(['\nMonitoring/original checkfit peak:\n' ...
            num2str(round(100*(monit_peak-check_resp_peak))/100) ' dB']);
        
        set(waveLine,'YData',monit)
        set(fftLine,'YData',20*log10(mon_spec/20e-6))
        
        inboxdisp=sprintf([mchpeak '\n' mck]);
        set(tekst,'String',inboxdisp)

        if isempty(A_sum), A_sum=0; end
        if isempty(B_sum), B_sum=0; end
        AB_corr=corr(A_sum,B_sum);
        drawnow
    end
    
    %% Record the response
    if k>measurements-10
        if round(i/2) == i/2
            stimLeft=measurements-k;
            IMPULSES=-repmat(impulses,1,stimLeft);
        else stimLeft=measurements-k;
            IMPULSES = repmat(impulses,1,stimLeft);
        end
    else
        if round(i/2) == i/2
            IMPULSES=-repmat(impulses,1,10);
        else IMPULSES=repmat(impulses,1,10);
        end
    end
    
    recDuration = length(IMPULSES)+delay_samples;
    pageNumber = playrec('playrec',IMPULSES', playChanList, recDuration,...
        recChanList);
    playrec('block',pageNumber);
    [MicIN] = playrec('getRec', pageNumber);
    
%     [correlation lags] = xcorr(MicIN,IMPULSES);
%     [CorMax Idx] = max(correlation);
%     delay_Samples = lags(Idx)

% filter the response between f1 and f2 
    MicIN=filter(coefb,coefa,MicIN);
% cut out the sound card delay
    MicIN=MicIN(floor(delay_samples)+1:end);
    

    %% Average 4 responses (linear cancellation) & store

    if k>measurements-10
       nm=measurements-k;
    else nm=10;
    end
    
    for m=1:nm
        response_1 = MicIN(1+(m-1)*length(MicIN)/nm:length(MicIN)*(m-3/4)/nm);
        response_1(1:round(time_rejection*fs))=0;
        response_2 = MicIN(1+length(MicIN)*(m-3/4)/nm:length(MicIN)*(m-1/2)/nm);
        response_2(1:round(time_rejection*fs))=0;
        response_3 = MicIN(1+length(MicIN)*(m-1/2)/nm:length(MicIN)*(m-1/4)/nm);
        response_3(1:round(time_rejection*fs))=0;
        response_4 = MicIN(1+length(MicIN)*(m-1/4)/nm:length(MicIN)*m/nm);
        response_4(1:round(time_rejection*fs))=0;
        
        if round(i/2) == i/2
            av_response=-(response_1+response_2+response_3+response_4)/4;
        else
            av_response=(response_1+response_2+response_3+response_4)/4;
        end
    
        av_response=av_response/(scIN*mic_sens); %[Pa]
        av_RMS=20*log10(sqrt(mean(av_response.^2))/20e-6); %[dB SPL]
        
        if av_RMS > threshold
            reject_count=reject_count+1;
        else k=k+1;
            stimulus(1:length(impulses),k)=MicIN((m-1)*length(impulses)+1:m*length(impulses));
        % store in 2 buffers, A and B
            if k/2==round(k/2); 
                B(:,k/2)=av_response;
            else
                A(:,(k+1)/2)=av_response;
            end
    
        end
    end
    
    B_sum=sum(B,2);
    set(BbLine,'YData',B_sum)
    A_sum=sum(A,2);
    set(AaLine,'YData',A_sum)
    
    Am = mean(A,2);
    Bm = mean(B,2);
    AB = (Am + Bm)./2;
    Noi = (Am - Bm)./2;
    WavRep = corr(Am,Bm);
    nfft_=size(Am,1);     
    AB_RMS = 20*log10(sqrt(mean(AB).^2)/20e-6);
    Noi_RMS = 20*log10(sqrt(mean(Noi).^2)/20e-6);
    AB_FreqAmp = 2*abs(fft(AB))/nfft_;
    AB_FreqAmp = 20*log10(AB_FreqAmp(1:nfft_/2)/20e-6);
    Noi_FreqAmp = 2*abs(fft(Noi))/nfft_;
    Noi_FreqAmp = 20*log10(Noi_FreqAmp(1:nfft_/2)/20e-6);
    
    set(TESpec,'ydata',AB_FreqAmp,'visible','on')
    set(NoiSpec,'ydata',Noi_FreqAmp,'visible','on')
    
    drawnow
    
    % calculate the percentage of rejected measurements
    reject_per=100*reject_count/i;
    rejected=['Rejected responses: ' num2str(reject_count)];
    AB_correlation=['AB waveforms correlation: ' num2str(round(100*AB_corr)/100) ...
      ' - '  num2str(WavRep)];
    LevText = sprintf('TEOAE level: %3.1f \nNoise level: %3.1f',...
       roundn(AB_RMS,-1),roundn(Noi_RMS,-1));
    
    set(rejtekst,'String',rejected)
    set(corrtekst,'String',AB_correlation)
    set(leveltekst, 'String',LevText)
    
    if reject_per > 50
        set(rejtekst,'ForegroundColor','r')
    end
    drawnow
    
end
% disp('End of measurement')
% pause(1)
% teoae_save1;

save_measurement(ID,Ear,'TE',...
                 {A,B,stimulus,Ear,PatientID},...
                 {'A','B','Stimulus','Ear','PatientID'},...
                 {fs,ydb,measurements,time_rejection,...
                  threshold,[f1*fs/2 f2*fs/2]},...
                 {'Fs','TEdBPeak','NumofMeasurements',...
                  'TimeRejection','NoiseThreshold',...
                  'BandPassFilter'});

maingui_new

% if(boolean_teoae_Proc)
%     close(niceFigure)
%     teoae_processing;
%     if(boolean_teaoe_Plot)
%         teoae_plotting;
%     else teoae_save;
%     end
% end
end
