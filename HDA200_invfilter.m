function Hd = HDA200_invfilter(inputsignal,cutoff)
% outputsignal = HDA200_invfilter(inputsignal,cutoff)
%Headphone equalization filter for the Sennheiser HDA-200 audiometric 
%headphones. The filter is based on the minimum phase inverse of the 
%averaged measured transferfuctions across 28 human subjects 
%(Hammershoi 1999). 
%All measured transferfuctions are normalized at 1kHz and averaged. 
%Then the minimum phase aproximation is obtained and inverted in the
%frequency domain. The amplitude of the frequency domain inverse responce
%is used to derrive a FIR filter using the frequency sampling method. The
%filter is limited at the low frequencies at the point where the inverse
%response starts to increase >0 dB (around 66 Hz). For these frequencies the
%filter is made to have a flat response. The high frequency limit is left 
%as a desig criteria. From the frequency specified by the "cutoff" input 
%the filter will act as a low-pass filter reaching an attenuation of 60 dB 
%at the niquist frequency. The filter is designed to work for signals with
%a 48000Hz sampling frequency.
%This fuctions need the file ptf5s09.mat
%Input arguments:
%   inputsignal:    time signal in Pa with a sampling frequency of 48kHz
%                   or
%                   the number of points for the calculation of the 
%                   frequency responce (real positive scalar). for the
%                   entire unit circle (0 to fs)
%   cutoff:         high-cutoff frequency in Hz, from this frequency and 
%                   above the filter acts as a low-pass filter.                
%Output Arguments:
%   outputsignal:   filtered version of the input time signal (Pa/samples)
%                   or 
%                   Frequency responce of filter in dB, of the length of 
%                   the scalar used as imput and plot of the measured 
%                   responces, the average, the ideal inverse and the 
%                   implemented filter. 
%
% Author:   Rodrigo Ordoñez
%           Acoustics, Department of Electronic Systems
%           Aalborg University
% Date:     28-08-2008
% Version:  1
%**************************************************************************

%%%%%%%%%%%%variables%%%%%%%%%%%%%%%%%%%% 
FiltOrder = 400;
fs = 48000;
if max(length(inputsignal))>50
    N = fs;%4096;%2^16; %fft sise
    freq = 0:fs/N:fs-fs/N';
else
    N = inputsignal;
    freq = 0:fs/inputsignal:fs-fs/inputsignal';
end
%%%%%%%%Hammersh�i 1999%%%%%%%%%
load ER_10C_inCoupler_ImpulseResponse

ptf = double(h(:,1));
Ptf = fft(ptf,N);
%Normalization at 1kHz;
Ptf = Ptf./Ptf(find(fix(freq)<=1000,1,'last'));
ptf = ifft(Ptf);
%%%minphase
[ptf_m Ptf_m] = minphase(ptf,N);
%%invert freq resp
Ptf_i = 1./Ptf_m;
PTF_i = abs(Ptf_i);

Lowfreqlim = find(20*log10(PTF_i)<0,1); % low cut where inverse function 
                                        %gives a low freq boost.

High1 = find(fix(freq)<=cutoff,1,'last'); 
High2= find(abs(20*log10(PTF_i(High1:end)))<=0.1,1); 
Highfreqlim = High1 + High2;
%high cut-off for filter... above 
                                      %this freq. the filter will act as a 
                                      %low pass filter
Niquist = find(fix(freq)<=24000,1,'last');
ramp = linspace(PTF_i(Highfreqlim+1),10^(-60/20),Niquist-Highfreqlim);

NormFreqs = freq(1:Niquist)./fs.*2;
Amplitudes = [ones(Lowfreqlim-1,1);PTF_i(Lowfreqlim:Highfreqlim); ones(Niquist-Highfreqlim,1)];

h = fdesign.arbmag('N,F,A', FiltOrder, NormFreqs, Amplitudes);

Hd = design(h, 'freqsamp');

%low pass filter

Fpass = 8000;       % Passband Frequency
Fstop = 16000;       % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 40;          % Stopband Attenuation (dB)
match = 'stopband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h1  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, fs);
Hd1 = design(h1, 'butter', 'MatchExactly', match);

if max(length(inputsignal)) > 50
 
    headphone_eq = filter(Hd,inputsignal);
    outputsignal = filter(Hd1,headphone_eq);
else
    headphone_eq = freqz(Hd,inputsignal,'whole',fs);
    lowpass_filt = freqz(Hd1,inputsignal,'whole',fs);
    outputsignal = 20*log10(abs(headphone_eq.*lowpass_filt));
    figure(1)
    % subplot(2,1,1)
    hold on
    set(gca,'box','on','linewidth',2)
    grid on
    axis([20 fs/2 -40 40])
    ax=axis;
    nr_of_dek=log10(ax(2)/ax(1));
    dB_per_dek=(ax(4)-ax(3))/nr_of_dek;
    whished_ratio=25;  %dB_per_dek.
    set(gca,'PlotBoxAspectRatio',[whished_ratio/dB_per_dek 1 1])
    set(gca,'XTickLabel',{'100' '1000' '10000'})
    set(gca,'XScale','log')
    set(gca,'Box','on')
    %plot(freq,20*log10(abs(PTF_s)),'g','LineWidth',0.000001)
    plot(freq,20*log10(abs(Ptf)),'r','LineWidth',2);
    plot(freq,20*log10(abs(Ptf_i)),'m','linewidth',2)
    plot(freq,outputsignal,'b','linewidth',2)
    plot(freq,20*log10(abs(headphone_eq)),'c','linewidth',2)
    xlabel('Frequency [Hz]','Fontsize',12)
    text(ax(1),ax(4)+5,'Gain [dB]','FontSize',12)
end
