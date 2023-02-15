function [Dp_lev Pri_lev Freqs TS TS_f fs T] = dpcalc(filename)
% % [Dp_lev Dp_noise Dp_std] = dpcalc(filename) 
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

load(filename);
T = Config.Duration;
Ndp = Config.NumofPoints;
Nreps = Config.NumofAvg;
fs = Config.Fs;
f1 = Data.F1Frequency;
f2 = Data.F2Frequency;
ratio = Config.Ratio;
fdp = Data.DPFrequency;
TS = Data.TimeSignal;%(0.15*fs:(T-0.15)*fs-1,:);
T = length(TS)/fs;

freqs = 0:1/T:fs-1;
%frequency domain:
TS_f = 2.*abs(fft(TS))./length(TS);
NFFT = length(TS);
% erbcoefs = MakeERBFilters(fs,fdp);
% for i = 1:Ndp
%     for j = 0:Nreps-1
%         TS_ERBs= ERBFilterBank(TS(:,i+j*16),erbcoefs);
%         TS_filt(:,i,j+1) = TS_ERBs(i,:)';
%     end
% end

%only include 1 ERB centered at the dp frequency.
%noise: energy average (fft amplitude in Pa) of the ERB 
%without the dp frequency
%dp level: amplitud of dp frequency
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
%     noiseh = TS_f(freqs>fdp(i)&freqs<=fu(i),i+j*16);
%     noisel = TS_f(freqs>=fl(i)&freqs<fdp(i),i+j*16);
%     noise(i,j+1) = mean([noisel; noiseh]);
%     dplev(i,j+1) = TS_f(freqs==fdp(i),i+j*16);
    
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

%%
% fig1 = figure('Position', [0 0 1400 900]);
% ax1 = axes('Parent',fig1);
% for k = 1:16
% plot(ax1,freqs(1:length(TS)/2),20*log10(TS_f(1:length(TS)/2,k:16:end)./20e-6))
% hold on 
% plot(ax1,F1(k,1),20*log10(f1_lev(k,:)./20e-6),'*','markersize',8)
% plot(ax1,F2(k,1),20*log10(f2_lev(k,:)./20e-6),'s','markersize',8)
% plot(ax1,FDP(k,1),20*log10(dplev(k,:)./20e-6),'o','markersize',8)
% axis(ax1,[800 5500 -30 70]);
% %set(gca,'XTick',f1l(k):50:f1u(k))
% hold off
% 
% M(k) = getframe;
% pause
% end
%movie(M,5,3)

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