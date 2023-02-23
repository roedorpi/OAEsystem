function [varargout] = dpcalc(filename,varargin)
% % [handles] = dpcalc(filename,handles) 
% Calculates mean level, and mean background noise with std in SPL using
% the recorded time signal from DPOAE measurements.
% Outputs: 
%Dp_lev: nX4 matrix Dim1: average DP level across repetitions; Dim2:
%DP standard deviation; Dim3: Average Noise level across 
%repetitions; Dim4: Noise standard deviation.
%
%Pri_lev:nX4matrix Dim1: Average level of F1 across repetitions; Dim2: F1 
%level standard deviation; Dim3 Average level of F2 across repetitions; 
%Dim4: F2 level standard deviation;
%
%Freqs: nX3 matrix Primary and dp frequencies (f1, f2, fdp) 
%

Config = filename.fio.Config;
Data = filename.fio.Data;
Ndp = Config.NumofPoints;
Nreps = Config.NumofAvg;
fs = Config.Fs;
f1 = Data.F1Frequency;
f2 = Data.F2Frequency;
ratio = Config.Ratio;
fdp = Data.DPFrequency;
TS = Data.TimeSignal;%(0.15*fs:(T-0.15)*fs-1,:);
T = length(TS)/fs;
%%bandpass filter 600 and 6000Hz
n = 6; Wn = [600 6000]*2/fs;
ftype = 'bandpass';
[b,a] = butter(n,Wn,ftype);
TS = filtfilt(b,a,TS);

% %%window time signal to avoid transients at start and finish
% Win = tukeywin(length(TS),0.2);
% for i = 1:size(TS,2);
%     TS(:,i) = TS(:,i).*Win;
% end
freqs = 0:1/T:fs-1/T;
%frequency domain:
TS_f = 2.*abs(fft(TS))./length(TS);

%only include 1 ERB centered at the dp frequency.
%noise: energy average (fft amplitude in Pa) of the ERB 
%without the dp frequency

ERB = 24.7*(4.37.*fdp./1000+1);
fu = (ERB + sqrt(ERB.^2+4.*fdp.^2))./2;
fl = (-ERB + sqrt(ERB.^2+4.*fdp.^2))./2;
%search for f1 in 1ERB centered at f1
ERB = 24.7*(4.37.*f1./1000+1);
f1u = (ERB + sqrt(ERB.^2+4.*f1.^2))./2;
f1l = (-ERB + sqrt(ERB.^2+4.*f1.^2))./2;
%search for f2 in 1ERB centered at f2
ERB = 24.7*(4.37.*f2./1000+1);
f2u = (ERB + sqrt(ERB.^2+4.*f2.^2))./2;
f2l = (-ERB + sqrt(ERB.^2+4.*f2.^2))./2;
for i = 1:Ndp
    for j = 0:Nreps-1
    
    tmpa = TS_f(freqs>f1l(i)&freqs<f1u(i),i+j*16);
    tmpb = freqs(freqs>f1l(i)&freqs<f1u(i));
    [f1_lev(i,j+1) indxa] = max(tmpa);
    F1(i,j+1) = tmpb(indxa);
    
    tmpc = TS_f(freqs>f2l(i)&freqs<f2u(i),i+j*16);
    tmpd = freqs(freqs>f2l(i)&freqs<f2u(i));
    [f2_lev(i,j+1) indxc] = max(tmpc);
    F2(i,j+1) = tmpd(indxc);
    
    FDP(i,j+1) = 2*F1(i,j+1)-F2(i,j+1);
    
    tmpe = TS_f(freqs>fl(i)&freqs<fu(i),i+j*16);
    tmpf = freqs(freqs>fl(i)&freqs<fu(i));
    dplev(i,j+1) = TS_f(roundn(freqs,-6)==roundn(fdp(i),-6),i+j*16);
    
    noiseh = TS_f(freqs>fdp(i)&freqs<=fu(i),i+j*16);
    noisel = TS_f(freqs>=fl(i)&freqs<fdp(i),i+j*16);
    noise(i,j+1) = mean([noisel; noiseh]);
%     ll = find(freqs<fdp(1),50,'last');
%     uu = find(freqs>fdp(1),50,'first');
%     noiseh = TS_f(uu(1:10:end),i+j*16);
%     noisel = TS_f(ll(1:10:end),i+j*16);
%     noise(i,j+1) = mean([noisel; noiseh]);
    
    end
end
%SPL averages across repetitions:
Dp_lev(:,1) = mean(20*log10(dplev./20e-6),2);
Dp_lev(:,3) = mean(20*log10(noise./20e-6),2);
Pri_lev(:,1) = mean(20*log10(f1_lev./20e-6),2);
Pri_lev(:,3) = mean(20*log10(f2_lev./20e-6),2);
%standard deviation of mean SPL
Dp_lev(:,2) = std(20*log10(dplev./20e-6),[],2);
Dp_lev(:,4) = std(20*log10(noise./20e-6),[],2);
Pri_lev(:,2) = std(20*log10(f1_lev./20e-6),[],2);
Pri_lev(:,4) = std(20*log10(f2_lev./20e-6),[],2);
%frequencies:
Freqs = cat(2,F1(:,1),F2(:,1),fdp');
% %set variables into handles
varargin{1}.Dp_lev = Dp_lev;
varargin{1}.Pri_lev = Pri_lev;
varargin{1}.Freqs = Freqs;
% varargin{1}.TS = TS;
% varargin{1}.TS_f = TS_f;
% varargin{1}.Leves(:,:,1) = f1_lev;
% varargin{1}.Leves(:,:,2) = f2_lev;
varargin{1}.Levels(:,:) = dplev;
varargin{1}.Noise(:,:) = noise;
varargin{1}.T = T;
varargin{1}.fs = fs;

varargout{1} = varargin{1};
%%

% fig1 = figure('Position', [0 0 800 600]);
% ax1 = axes('Parent',fig1);
% set(ax1,'fontsize',14, 'box','on')
% xlabel(ax1,'Frequency [Hz]','fontsize',14);
% ylabel(ax1,'SPL [dB re. 20\mu Pa]','fontsize',14);
% ax = [500 5500 -30 75];
% axis(ax1,ax);
% 
% for k = 1:16
%     set(ax1,'nextplot','add')
%     plot(ax1,freqs(1:length(TS)/2),20*log10(TS_f(1:length(TS)/2,k)./20e-6))
%     plot(ax1,F1(k,1),20*log10(f1_lev(k,1)./20e-6),'gd',...
%         'markersize',8, 'markerfacecolor','g')
%     plot(ax1,F2(k,1),20*log10(f2_lev(k,1)./20e-6),'ms',...
%         'markersize',8, 'markerfacecolor','m')
%     plot(ax1,FDP(k,1),20*log10(dplev(k,1)./20e-6),'ro',...
%         'markersize',8, 'markerfacecolor','r')
%     M(k) = getframe(fig1);
%     cla(ax1)
% end
% M = repmat(M,1,5);
% movie2avi(M,'dpmeanM05.avi','compression','none','fps',4,'quality',100);
% close(fig1)
% movie(M,1,5)

% 
% hold on
% axis([700 3500 -20 25])
% box on
% grid on
% area(fdp,Dp_noise,-20,'facecolor',[0.7 0.7 0.7])
% plot(fdp,Dp_lev,'bo-','linewidth',2,'markersize',8,'markerfacecolor','b')
% errorbar(fdp,Dp_lev,Dp_std,'b-')


% ERB plot from AuditoryToolbox
%   y = ERBFilterBank([1 zeros(1,NFFT-1)], fcoefs);
% 	resp = 20*log10(abs(fft(y')));
% 	freqScale = (0:NFFT-1)/NFFT*fs;
% 	semilogx(freqScale(1:NFFT/2-1),resp(1:NFFT/2-1,:));
% 	axis([100 16000 -60 0])
% 	xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');