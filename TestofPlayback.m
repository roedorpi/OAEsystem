%% 
clear
clc
close all

%%
Main = figure(1);
ax = axes('Parent',Main,'Tag','AxTime','box','on');
ax1 = axes('Parent',Main,'Tag','AxTime1','box','on');
btn = uicontrol('Parent',Main,'Style','pushbutton','String','CheckFit','Callback',@runcheck);
ax.Units ='normalized';

btn.Units = 'normalized';
ax.Position = [0.15 0.57 0.8 0.4];
ax1.Position = [0.15 0.15 0.8 0.4];
btn.Position = [0.05 0.05 0.1 0.05];

%%
AuIO = audioPlayerRecorder;
b = getAudioDevices(AuIO);
AuIO.Device = b{2};
AuIO.SampleRate = 48000;
AuIO.PlayerChannelMapping = [1];
AuIO.RecorderChannelMapping = [1];
% Check the buffer size properties on the ASIO settings to set the same
% here.
fs = AuIO.SampleRate;
AuIO.SupportVariableSize = true;
AuIO.BufferSize = 512;

%%
%dvs = playrec('getDevices');
playrec('reset');
playrec('init',fs,7,7,1,1);

%%


scDAC = 1;
check_w=80e-6;
check_level = 80;
check_l = 0.1;

check_Pa=20*10^((check_level/20)-6);
check_digi=check_Pa/scDAC;

spacing = 20e-3;
yPa=20*10^((check_level/20)-6);
imp_level=yPa/scDAC; %[digital]

check_impulse=ones(1,round(check_w*fs));
check_impulse=[zeros(1,100) check_impulse zeros(1,round((check_l-check_w)*fs)-100)]; %impulse length 0.1 s

% compensating for the loudspeaker transfer function
% check_impulse = lsinvopt10kHz(check_impulse')';

% low pass filtering
load('lowfilter20db.mat')
check_impulse=filter(hd,check_impulse);

% adjust the level after filtering
check_impulse_peak = max(abs(check_impulse));
correction = check_digi/check_impulse_peak;
check_impulse=check_impulse.*correction;


impulse1 = zeros(1,spacing*fs);
impulse1(1:round(check_w*fs)) = 1/3;

impulse3 = zeros(1,spacing*fs);
impulse3(1:round(check_w*fs)) = -1;

impulses = [impulse1 impulse1 impulse1 impulse3];

%compensating for the loudspeaker transfer function
%impulses = lsinvopt10kHz(impulses')';

% low pass filtering
impulses=filter(hd,impulses);

% adjust the level after filtering
impulses=imp_level*impulses/max(abs(impulses)); %0.8515;

impulses = [zeros(1,round((check_l - length(impulses)/fs)*fs)) impulses ]; 
impulses = [repmat(impulses,1,10) zeros(1,length(impulses))];

%%
delete(ax.Children)
ax.NextPlot = "add";
x = playRecSig(AuIO,zeros(length(check_impulse),1));
t = 1/fs:1/fs:length(x)/fs;

pl1 = plot(ax,t,x(:,1));
t_ = 1/fs:1/fs:length(check_impulse)/fs;
pl2 = plot(ax,t_,check_impulse);

delete(ax1.Children)
x1 = playRecSig(AuIO,zeros(length(impulses),1));
t1 = 1/fs:1/fs:length(x1)/fs;
pl3 = plot(ax1,t1,x1);
hold on
t_ = 1/fs:1/fs:length(impulses)/fs;
pl4 = plot(ax1,t_,impulses);


%%
delete(ax.Children)
ax.NextPlot = "add";
out = playrec('playrec',zeros(length(check_impulse),1),1,length(check_impulse)+20e-3*fs,1);
while(playrec('isFinished', out) == 0);end
x = playrec('getRec',out);
t = 1/fs:1/fs:length(x)/fs;
pl1 = plot(ax,t,x(:,1));
t_ = 1/fs:1/fs:length(check_impulse)/fs;
pl2 = plot(ax,t_,check_impulse);


delete(ax1.Children)
ax1.NextPlot = "add";
out = playrec('playrec',zeros(length(impulses),1),1,2*length(impulses),1);
while(playrec('isFinished', out) == 0);end
x = playrec('getRec',out);
t = 1/fs:1/fs:length(x)/fs;
pl3 = plot(ax1,t,x(:,1));
t_ = 1/fs:1/fs:length(impulses)/fs;
pl4 = plot(ax1,t_,impulses);

%%
T = 0;
tStart = tic;
while T < 30
    T = toc(tStart);
    
%     out = playrec('playrec',impulses',1,2*length(impulses),1);
%     while(playrec('isFinished', out) == 0);end
%     x_te = playrec('getRec',out);

    x_te = playRecSig(AuIO,impulses');
    x_cf = playRecSig(AuIO,check_impulse');
    %x_ = playRecSig(AuIO,impulses');
    
    %x = mean(cat(3,x,x_),3);
    pl1.YData = -x_cf(:,1);
    pl3.YData = -x_te(:,1);
    drawnow
end

