close all 
clear all
clc



load ER_10C_inCoupler_ImpulseResponse

h_ = h(2000:2800,:);
t = 0:1/fs:length(h_)/fs -1/fs;

winL=0.003; % 40ms length
tt=0:1/fs:winL-1/fs;
win=hanning(length(tt));
Window=[win(1:length(win)/2); ones(length(t)-length(win),1); ...
win((length(win)/2)+1:end)];
            
y = h_.*repmat(Window,1,4);




nfft = 4800;

Y = 20*log10(abs(fft(y,nfft)));
Y(:,1) = Y(:,1)-Y(101,1);
Y(:,2) = Y(:,2)-Y(101,2);
Y(:,3) = Y(:,3)-Y(101,3);
Y(:,4) = Y(:,4)-Y(101,4);
F = 0:fs/nfft:fs-fs/nfft;




fig = figure(1);

ax1 = axes('Parent',fig,'position',[0.1 0.6 0.8 0.35],'box','on');
ax2 = axes('Parent',fig,'position',[0.1 0.1 0.8 0.35],'box','on');
plot(ax1,t,y)
plot(ax2,F(1:round(end/2)),Y(1:round(end/2),:));

xlabel(ax1, 'Time [s]');
ylabel(ax1, 'Amplitude');
set(ax2,'xscale','log','xlim',[20 22000],'ylim',[-30 30])

xlabel(ax2, 'Frequency [Hz]');
ylabel(ax2, 'Amplitude [dB]');


