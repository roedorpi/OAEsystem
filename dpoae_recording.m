function dpoae_recording()
global type_range ratio f1min f1max num_points level_f1 level_f2 num_av ...
    duration type_dpoae fs...
    temp freq_f1 freq_f2 noise_av dp_level dp_freq ID Ear...
    mic_sens scIN scLOUT scROUT speakerL speakerR...
    delay_samples


    

%% Instruction

%% 
%This function send the stimuli to the loudspeaker and record the ear
%response.

% This function requires a running playrec installation, follow installation
%instructions at http://www.playrec.co.uk/

% This function plot the spectogramme of what is recorded by the
% microphone, calculate the dpoae for a range of frequencies, estimate the
% noise, plot and save them in matrices call dpoae and noise_est (each
% column is for a new subject).

%Parameters:
%type_range can be 'logari' if frequency scale for dpoae wanted is
%   logarithmic or 'linear' if the frequency scale wanted is linear
%ratio is the ratio between f1 and f2, f2 = ratio*f1
%f1min is the minimum frequency wanted for f1
%f1max is the maximum frequency wanted for f1
%num_points is the number of points wanted between f1min and f1max
%level_f1 and level_f2 are the levels of the stimulis to send
%num_av is the number of repetitions.
%duration is the duration of each pair of stimulis.
%type_dpoae is the type of dpoae wanted to be measure for example 2*f2-f1
%subject_number is the number of the subject

%% Initialisation
%Frequency range for the spectogramme
f = 0:1/duration:(fs/2);
% Length of the fft
L = 2*length(f);
%Frequency range for the dpoae
 for k = 0:num_points-1
    freq(1,k+1) = k;
 end

% Define the different frequencies to send depending of the parameters type_range, f1min, f1max and ratio     
if strcmp(type_range,'linear') == 1
    sample_f1 = (f1max-f1min)/(num_points-1);
    f1 = f1min:sample_f1:f1max;
elseif strcmp(type_range,'logari')
    for k = 0:num_points-1
        f1(1,k+1) = (f1max/f1min)^(k/(num_points-1))*f1min;    
    end
else
    display('Error in DPOAE run. Unknown f1 spacing argument.')
    quit;
end
f2 = ratio.*f1;
freq_f1=f1; freq_f2=f2;

% Definition of f3, frequency of the dpoae depending of type_dpoae

freqs = 0:1/duration:fs-1;

if strcmp(type_dpoae,'2f1-f2') == 1
    f3 = 2.*f1-f2;
    resolution = 1/duration;
    fd = round((f3)/resolution);
    f3 = freqs(fd);
    f2 = 2*f1-f3;%*resolution;
%     f3 = 2.*f1-f2;
    f2max = max(f2);
    f2min = min(f2);
    f3max = 2*f1max-f2max;
    f3min = 2*f1min-f2min;
elseif strcmp(type_dpoae,'3f1-2f2') == 1
    f3 = 3.*f1-2*f2;
    resolution = 1/duration;
    fd =round((f3)/resolution);
    f2 = (3*f1-fd*resolution)/2;
    f3 = 3*f1-2*f2;
    f2max = max(f2);
    f2min = min(f2);
    f3max = 3*f1max-2*f2max;
    f3min = 3*f1min-2*f2min;
else
    display('Error in DPOAE run. Unknown DP frequency.')
    quit;
end

dp_freq=f3;

% Amplitude of the signals to send
amplitude1 = zeros(1,num_points);
amplitude2 = zeros(1,num_points);
for i=1:num_points
    f1_level = DPSPLCorrection(level_f1,f1(1,i),1);
    f2_level = DPSPLCorrection(level_f2,f2(1,i),2);

    amplitude1(1,i) = sqrt(2)*10^(f1_level/20)*20*10^(-6)/(speakerR*scROUT);
    amplitude2(1,i) = sqrt(2)*10^(f2_level/20)*20*10^(-6)/(speakerL*scLOUT);
    if  amplitude1(1,i)>1 
        amplitude1(1,i) = 1;
        display('Warning! F1 level is to high')
    end    
    if  amplitude2(1,i)>1 
        amplitude2(1,i) = 1;
        display('Warning! F2 level is to high')
    end
end

% Vector time for the signals to send
t = 0:1/fs:duration-1/fs;

stimuli1 = zeros(length(t)+6000,num_points);
stimuli2 = zeros(length(t)+6000,num_points);

% time window for the stimuli
winL=0.1; % 100ms length
tt=0:1/fs:winL-1/fs;
win=sin(2*pi*tt*(1/winL/2));
window=[win(1:length(win)/2) ones(1,length(t)-length(win)) win((length(win)/2)+1:end)];

for i=1:1:num_points
stimuli1(1:length(t),i) = amplitude1(1,i).*sin(2*pi.*f1(1,i).*t).*window; 
stimuli2(1:length(t),i) = amplitude2(1,i).*sin(2*pi.*f2(1,i).*t).*window;
% if i==num_points
%     uisave('stimuli1')
% end
end

%Initialisation of the vector which are going to contain the average values
dpo_av = zeros(1,num_points);
noise_av = zeros(1,num_points);

% Initialisation of the figure for the spectogramme
fftFigure = figure;
dcm_obj = datacursormode(fftFigure);
set(fftFigure,'ToolBar','none','MenuBar','none');
scrsize = get(0,'ScreenSize');
set(fftFigure,'OuterPosition',[0 30 scrsize(3) scrsize(4)-30])
fftAxes = axes('parent', fftFigure, 'xlimmode', 'manual', 'ylimmode', 'manual', 'xscale', 'log', 'yscale', 'linear', 'xlim', [10 fs/2], 'ylim', [-20, 80]);
fftLine = line('XData',f,'YData', f);
set(fftAxes,'Position',[0.05 0.6 0.9 0.35]) 
xlabel(fftAxes,'Frequency (Hz)')
ylabel(fftAxes,'dB SPL')
title(fftAxes,'FFT of the recorded signal')
drawnow;

%Initialisation of the figure for the dpoae if type_range is logari
if(strcmp(type_range,'logari'))
dpoaeAxes = axes('parent', fftFigure, 'xlimmode', 'manual', 'ylimmode', 'manual', 'xscale', 'log', 'yscale', 'linear', 'xlim', [round(f3min) round(f3max)], 'ylim', [-20, 20]);
grid on;
dpoaeLine = line('XData', f3,'YData', f3,'Marker','o');
dpoaeLine2 = line('XData', f3,'YData', f3,'Marker','o');
set(dpoaeAxes,'Position',[0.05 0.1 0.9 0.4]) 
set(dpoaeAxes,'Xtick',round(f3));
xlabel(dpoaeAxes,'Dpoae Frequency (Hz)')
ylabel(dpoaeAxes,'dB SPL')
title(dpoaeAxes,'DPOAE and Noise Estimation')
legend('Noise Estimation','DPOAE')
drawnow;
end
%Initialisation of the figure for the dpoae if type_range is linear
if(strcmp(type_range,'linear'))
dpoaeAxes = axes('parent', fftFigure, 'xlimmode', 'manual', 'ylimmode', 'manual', 'xscale', 'linear', 'yscale', 'linear', 'xlim', [round(f3min)-10 round(f3max)], 'ylim', [-20, 20]);
grid on;
dpoaeLine = line('XData', f3,'YData', f3,'Marker','o');
dpoaeLine2 = line('XData', f3,'YData', f3,'Marker','o');
set(dpoaeAxes,'Position',[0.05 0.1 0.9 0.4])
set(dpoaeAxes,'Xtick',round(f3));
xlabel(dpoaeAxes,'Dpoae Frequency (Hz)')
ylabel(dpoaeAxes,'dB SPL')
title(dpoaeAxes,'DPOAE and Noise Estimation')
legend('Noise Estimation','DPOAE')
end
drawnow;

%% Processing
temp = zeros(length(stimuli1(:,num_points))-6000,num_points*num_av);
for i = 1:num_av; %Number of averages
    for j = 1:1:num_points %Number of point for each average
        out = playrec('playrec', [stimuli1(:,j) stimuli2(:,j)] ,[1 2], length(stimuli1(:,j))-6000+delay_samples, 2);
        playrec('block',out);
        record = playrec('getRec',out)/mic_sens/scIN;  
%         delay_samples
%         size(record)
%         size(temp)
        record = record(delay_samples+1:end);
        temp(:,(i-1)*num_points+j) = record;
        L = length(record);
        fftrecordtp=zeros(1,L);
        fftrecord=zeros(1,L/2+1);
        fftrecordtp(1,:)= 2*abs(fft(record(:,1),L))/L;
        fftrecord(1,:) = fftrecordtp(1,1:round(L/(2))+1);
        % Make sure that the point is the right one (the one with the most
        % energy)
        [val,fd1]=max(fftrecord(1,fd(1,j)-1:fd(1,j)+1));
        fd2 = fd(1,j)+fd1-2;
        % Noise estimation
        noise = fftrecord(1,[fd2-5:fd2-2,fd2+2:fd2+5]);
        noise_mean(i,j)= 20*log10(mean(noise)/20e-6);
        % Dpoae estimation
        dpo(i,j) = 20*log10(fftrecord(fd2)/20e-6);
        %Dpoae and noise averaging
        dpo_av(1,j) = (dpo_av(1,j)*i+dpo(i,j))/(i+1);
        noise_av(1,j) = (noise_av(1,j)*i+noise_mean(i,j))/(i+1);
        %Plot spectogramme
        set(fftLine, 'YData', 20*log10(fftrecord(1,:)/20e-6));
        drawnow;
        %Plot dpoae and noise estimation 
        set(dpoaeLine, 'YData', noise_av,'Color','red');
        set(dpoaeLine2, 'YData', dpo_av,'Color','blue');
        drawnow;
    end
end

dp_level=dpo;

% uisave({'noise_av','dp_level','freq_f1','freq_f2','dp_freq'})

set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','on')
set(dcm_obj,'Enable','on','UpdateFcn',{@myupdatefcn,f3,dpo_av,noise_av,fftFigure})
%close(fftFigure);

save_measurement(ID,Ear,'DP',...
                 {level_f1,level_f2,freq_f1,freq_f2,dp_level,dp_freq,noise_av,temp},...
                 {'F1Level','F2Level','F1Frequency','F2Frequency','DPLevel','DPFrequency',...
                 'NoiseLevel','TimeSignal'},...
                 {fs,type_range,ratio,f1min,f1max,num_points,num_av,duration},...
                 {'Fs','FrequencyRange','Ratio','F1Min','F1Max','NumofPoints',...
                 'NumofAvg','Duration'});

maingui_new

end