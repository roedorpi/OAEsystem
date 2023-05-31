function teoae_recording_new(fig)
h = getappdata(fig,'main');
% This function masures TEOAE using non-linear mode of stimuli. Has to be 
% run after checkfit.m, because of the check-fit monitoring included in it.
% Input parameters are checkfit response, delay (from checkfit function),
% input and output channels, stimuli level, number of measurements,
% rejection time, rejection threshold, band pass filtering cutoff
% frequencies.
% The function outputs two matrices, with time signals from all the
% measurements. The TEOAE processing funcion can be run afterwards to get
% more TEOAE parameters (mean response, reproducibility, etc.).

% %global A B check_response resp_spec check_impulse check_level check_level_adj...
%     delay_samples checkfit_monit fs playChanList recChanList ydb ...
%     measurements time_rejection threshold spacing f1 f2 niceFigure monit ...
%     tt ABAx AaLine BbLine scLOUT scROUT scIN speakerL speakerR mic_sens ...
%     chwindow stimulus

% global PatientID Ear ID

load('lowfilter20db.mat')

Ear  = h.MeasParamTab.Data{11,2};
%%%%%%%%%%%%%%%%%%%%%%% scaling factors %%%%%%%%%%%%%%%%%%%%%%%%
% scLOUT      = 3.3220;   % [Vrms/output digital rms] output 1 (left speaker)
% scROUT      = 3.4266;   % [Vrms/output digital rms] output 2 (right speaker)
% speak_sens  = 0.07963;  % [Pa/Vrms]
% mic_sens    = 0.5;      % [Vrms/Pa]
% scIN        = 0.433;    % [digital rms/Vrms]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose channels 
%if playChanList == 1; scOUT = scLOUT; else scOUT = scROUT; end
%if playChanList == 1; speak_sens = speakerL; else speak_sens = speakerR; end

for i = 1:4, cla(h.ax(i)), end

scDAC = h.Config.CheckFit.scDAC;
scADC = h.Config.CheckFit.scADC;
fs = h.AuIO.SampleRate;
%% Generate impulse stimuli
check_level = h.MeasInfoTab.Data{8,2};
check_w = h.Config.CheckFit.check_w;
check_l = h.Config.CheckFit.check_l;
ydb = h.MeasParamTab.Data{1,2};
time_rejection = h.MeasParamTab.Data{3,2};
noise_rejection = h.MeasParamTab.Data{2,2};
spacing = h.MeasParamTab.Data{4,2}/1000;
measurements = h.MeasParamTab.Data{5,2};
fLow = h.MeasParamTab.Data{6,2};
fHigh = h.MeasParamTab.Data{7,2};
checkfit_monint = h.MeasParamTab.Data{9,2};
LevDiff = ydb - check_level;
StimMaxPeakLev =  ydb + LevDiff;

yPa=20*10^((StimMaxPeakLev/20)-6);
imp_level=yPa/scDAC; %[digital]

    impulse1 = zeros(1,spacing*fs);
    impulse1(1:round(check_w*fs)) = 1/3;

    impulse3 = zeros(1,spacing*fs);
    impulse3(1:round(check_w*fs)) = -1;

impulses = [impulse1 impulse1 impulse1 impulse3];

%compensating for the loudspeaker transfer function
impulses = lsinvopt10kHz(impulses')';

% low pass filtering
impulses=filter(hd,impulses);

% adjust the level after filtering
impulses=imp_level*impulses/max(abs(impulses)); %0.8515;

%impulses = [zeros(1,round((check_l - length(impulses)/fs)*fs)) impulses ]; 

% stimulus=impulses;

%% Probe fit response

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
        check_response = h.Data.ProbeFitL(:,end);
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
        check_response = h.Data.ProbeFitR(:,end);
end

% draw TEOAE resp
tt=(0:1/fs:((length(impulses)/4)-1)/fs)*1000;

AaLine = line('parent',h.ax(1),'XData',tt,'YData',zeros(1,length(tt)),'Color','b');
BbLine = line('parent',h.ax(1),'XData',tt,'YData',zeros(1,length(tt)),'Color','g');
legend(h.ax(1),'Averaged waveform A','Averaged waveform B','Location','NorthWest')


ff = (0:((length(impulses)/4))/2-1)/((length(impulses)/4)).*fs;
TESpec = area(h.ax(2),ff,zeros(1,length(ff)),-60,'facecolor','g',...
    'visible','off');
set(h.ax(2),'nextplot','add')
NoiSpec = area(h.ax(2),ff,zeros(1,length(ff)),-60,'facecolor','r',...
    'visible','off');

legend(h.ax(2),{'TEOAE' 'Noise'},'location','northeast','fontsize',14)
legend(h.ax(2),'boxoff')
set(h.ax(2),'ylim',[-60 10], 'xlim', [0 6000],'xtick',1000:1000:6000,...
    'xticklabel',{'1' '2' '3' '4' '5' ''},'xgrid','on','ygrid','on',...
    'box','on');

% draw Chefit resp. 
waveLine = line('parent',h.ax(3),'XData',t_,'YData',zeros(length(t_),1));
fftLine = area(h.ax(4),f_,zeros(1,length(f_)),0);


drawnow

k=0;
reject_count=0;
AB_corr=1;
stimulus=NaN(length(impulses),measurements);

%filtering the response (coefficients)
fLow=2*fLow/fs; fHigh=2*fHigh/fs;
[coefb, coefa]=butter(4,[fLow fHigh]);
% 
A_sum=[]; B_sum=[];
if round(checkfit_monint/10) ~= checkfit_monint/10
     checkfit_monint = ceil(checkfit_monint/10)*10;
end

i = 0;
while k<measurements
   %% Checkfit monitoring
   i = i + 1;
   if round((i-1)/(checkfit_monint/10))==(i-1)/(checkfit_monint/10) && i~=1
        monit = playRecSig(h.AuIO,[h.Config.CheckFit.check_impulse' ...
            zeros(length(h.Config.CheckFit.check_impulse),1)],0);
        
        monit=hpf(monit,500,fs);
        
        win=hann(400);
        chwindow=[win(1:length(win)/2); ones(length(check_response)-length(win),1);...
        win(length(win)/2+1:length(win))];
        monit=monit(:,1).*chwindow;
        
        monit=monit/scADC; %[Pa]
        monit_peak=20*log10(max(abs(monit))/20e-6);
        
        mon_spec=abs(fft(monit',nfft_))/length(monit);
        mon_spec=2*mon_spec(1:length(mon_spec)/2);
        mon_check_corr=corr(monit,check_response);
        % Check fit correlation for moniotoring probe placement. 
        h.MeasInfoTab.Data{10,2} = roundn(mon_check_corr,-2);
        check_resp_peak=20*log10(max(abs(check_response))/20e-6);
        % Peak ratio of initial fit and latest measured 
        h.MeasInfoTab.Data{11,2} = roundn(monit_peak-check_resp_peak,-1);
        
        set(waveLine,'YData',monit)
        set(fftLine,'YData',[mon_spec(1) 20*log10(mon_spec(2:end)/20e-6)])
       
        if isempty(A_sum), A_sum=0; end
        if isempty(B_sum), B_sum=0; end
        AB_corr=corr(A_sum,B_sum);
        drawnow
    end
    
    %% Record the response
%     if k>measurements-10
%         if round(i/2) == i/2
%             stimLeft=measurements-k;
%             IMPULSES=-repmat(impulses,1,stimLeft);
%         else stimLeft=measurements-k;
%             IMPULSES = repmat(impulses,1,stimLeft);
%         end
%     else
        if round(i/2) == i/2
            IMPULSES=-repmat(impulses,1,10);
        else IMPULSES=repmat(impulses,1,10);
        end
%     end
    IMPULSES = [IMPULSES zeros(1,1.5*length(impulses))]; %h.Config.CheckFit.scDelay;
    MicIN = playRecSig(h.AuIO,[IMPULSES',zeros(length(IMPULSES),1)],h.MeasParamTab.Data{10,2});
    InDXm = find(abs(MicIN(:,1))> 0.8*(max(abs(MicIN(:,1)))),1);
    InDX = find(abs(IMPULSES)> 0.8*(max(abs(IMPULSES))),1);
    scDelay = InDXm-InDX;
    h.Config.CheckFit.scDelay = mean([h.Config.CheckFit.scDelay scDelay]);
    %scDelay(i) = length(MicIN)-length(IMPULSES)
    
    if scDelay > 0 

    % filter the response between f1 and f2 
        MicIN=filter(coefb,coefa,MicIN(:,1));
    % cut out the sound card delay
        MicIN=MicIN(scDelay:end);
        MeasLength = length(MicIN) - scDelay;
        if MeasLength > spacing*fs*4
            TailEnd = rem(MeasLength,spacing*fs*4);
            nm = floor(MeasLength/(spacing*fs*4));
            MicIN = MicIN(1:spacing*fs*4*nm);
            %% Average 4 responses (linear cancellation) & store
        
        %     if k>measurements-10
        %        nm=measurements-k;
        %     else nm=10;
        %     end
            
            Resps = reshape(MicIN,spacing*fs*4,nm);
            Responses = reshape(MicIN,spacing*fs,4,nm);
            Responses(1:round(time_rejection/1000*fs),:,:) = 0;
            if round(i/2) == i/2
                AverResp = -squeeze(mean(Responses,2));
            else
                AverResp = squeeze(mean(Responses,2));
            end
            for m = 1:size(AverResp,2)
            
                av_response=AverResp(:,m)/(scADC); %[Pa]
                av_RMS=20*log10(sqrt(mean(av_response.^2))/20e-6); %[dB SPL]
                
                if av_RMS > noise_rejection
                    reject_count=reject_count+1;
                else k=k+1;
                    stimulus(1:length(impulses),k)= Resps(:,m);
                % store in 2 buffers, A and B
                    if k/2==round(k/2) 
                        B(:,k/2)=av_response;
                    else
                        A(:,(k+1)/2)=av_response;
                    end
            
                end
            end
        end
    end

    B_sum=sum(B,2);
    
    A_sum=sum(A,2);
   
    
    Am = mean(A,2);
    Bm = mean(B,2);
    set(BbLine,'YData',B_sum)
    set(AaLine,'YData',A_sum)
    
    AB = (Am + Bm)./2;
    Noi = (Am - Bm)./2;
    WavRep = corr(Am,Bm);
    nfft__=size(Am,1);     
    AB_RMS = 20*log10(sqrt(mean(AB).^2)/20e-6);
    Noi_RMS = 20*log10(sqrt(mean(Noi).^2)/20e-6);
    AB_FreqAmp = 2*abs(fft(AB))/nfft__;
    AB_FreqAmp = 20*log10(AB_FreqAmp(1:nfft__/2)/20e-6);
    Noi_FreqAmp = 2*abs(fft(Noi))/nfft_;
    Noi_FreqAmp = 20*log10(Noi_FreqAmp(1:nfft__/2)/20e-6);
    
    set(TESpec,'ydata',AB_FreqAmp,'visible','on')
    set(NoiSpec,'ydata',Noi_FreqAmp,'visible','on')
    
    drawnow
    
    % calculate the percentage of rejected measurements
    reject_per=100*reject_count/i;

    %rejected=['Rejected responses: ' num2str(reject_count)];
    h.MeasInfoTab.Data{6,2} = WavRep;
    h.MeasInfoTab.Data{7,2} = round(100*AB_corr)/100;

    h.MeasInfoTab.Data{3,2} = roundn(AB_RMS,-1);
    h.MeasInfoTab.Data{4,2} = roundn(Noi_RMS,-1);
    h.MeasInfoTab.Data{5,2} = reject_per;
    
%     if reject_per > 50
%         set(rejtekst,'ForegroundColor','r')
%     end
    drawnow
    
end
h.AuIO.release;
% disp('End of measurement')
% pause(1)
% % teoae_save1;
h.Data.A = A;
h.Data.B = B;
h.Data.Stimulus = stimulus;
h.Data.PatientID = h.File.Name(1:end-4);
h.Data.Date = h.MeasInfoTab.Data{1,2};
h.Data.StartTime = h.MeasInfoTab.Data{2,2};
h.Data.EndTime = datestr(rem(now,1),'HH:MM:SS');
h.Data.Ear = h.MeasParamTab.Data{11,2};

h.Config.Fs = fs;
h.Config.TEdBPeak = h.MeasParamTab.Data{1,2};
h.Config.NumofMeasurements = size(h.Data.Stimulus,2);
h.Config.TimeRejection = h.MeasParamTab.Data{3,2};
h.Config.NoiseThreshold = h.MeasParamTab.Data{2,2};
h.Config.BandPassFilter = [h.MeasParamTab.Data{6,2} h.MeasParamTab.Data{7,2}];
h.Config.TEsignal = IMPULSES;
h.Config.CLS = h.MeasParamTab.Data{10,2};
TimeStamp = datestr(now,'mmm_dd_yyyy_HH_MM_SS');
h.FileIO.(TimeStamp) = struct('Data',h.Data,'Config',h.Config);


setappdata(fig,'main',h)

uiwait(msgbox('TEOAE data saved!'))


% save_measurement(ID,Ear,'TE',...
%                  {A,B,stimulus,Ear,PatientID},...
%                  {'A','B','Stimulus','Ear','PatientID'},...
%                  {fs,ydb,measurements,time_rejection,...
%                   threshold,[f1*fs/2 f2*fs/2]},...
%                  {'Fs','TEdBPeak','NumofMeasurements',...
%                   'TimeRejection','NoiseThreshold',...
%                   'BandPassFilter'});
% 
% maingui_new

% if(boolean_teoae_Proc)
%     close(niceFigure)
%     teoae_processing;
%     if(boolean_teaoe_Plot)
%         teoae_plotting;
%     else teoae_save;
%     end
% end
end
