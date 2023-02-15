function Hd = InvFilter(y,nfft,cutoff,FiltOrder,fs)
% Hd = InvFilter(InpulseResp,nfft,cutoff,FiltOrder,fs)
%Minimum phase inverse filter of transfer function measurement. 
%
% The measured transferfuction is normalized at 1kHz, then the minimum 
% phase aproximation is obtained and inverted in the frequency domain. The 
% amplitude of the frequency domain inverse is used to derrive a FIR filter 
% using the frequency sampling method. The filter is limited at the low 
% frequencies at the point where the inverse response starts to increase 
% >0 dB. For these frequencies the filter is made to have a flat response. 
% The high frequency limit is left as a desig criteria. From the frequency 
% specified by the "cutoff" input the filter will also have a flat response
% 
% This fuctions need the file minphase.mat
% Input arguments:
%   y:          time signal of the impulse response to invert.
%   nfft:       the number of points for the calculation of the frequency 
%               responce (real positive scalar), for the entire unit circle 
%               (0 to fs).
%   cutoff:     high-cutoff frequency in Hz, from this frequency and above 
%               the filter has a flat response.
%   FiltOrder:  Order of the desired inverse FIR filter
% 
% Output Arguments:
%   Hd:         Matlab Filter Structure to be used with the filter command. 
%
%                   filtered_signal = filter(Hd,unfiltered_signal)
%
%               Both input and output are time signals. 
%
% Author:   Rodrigo Ordoñez
%           Acoustics, Department of Electronic Systems
%           Aalborg University
% Date:     18/07/2012
%**************************************************************************

addpath /media/data/work/experiments/functions/
%% Frequency and time vectors
t = 1/fs:1/fs:nfft/fs;
freq = 0:fs/nfft:fs-fs/nfft';
%% Input impulse response
[Max I1d] = max(abs(y));
premax = 20;
resplength = 1024;
%y = y(I1d-premax:I1d+resplength-premax);

ptf = double(y(:));
Ptf = fft(ptf,nfft);
% Normalization at 1kHz;
Ptf = Ptf./Ptf(find(fix(freq)<=1000,1,'last'));
ptf = ifft(Ptf);
% minphase
%[ptf_m Ptf_m] = minphase(ptf,nfft);
% invert frequency response
Ptf_i = 1./Ptf;
PTF_i = abs(Ptf_i);

% low frequency cut-off where inverse function gives a low freq boost.
Lowfreqlim = find(20*log10(PTF_i)<0,1); 

% high frequency cut-off, above this freq. the filter has a flat response.
High1 = find(fix(freq)<=cutoff,1,'last'); 
High2 = find(abs(20*log10(PTF_i(High1:end)))<=0.5,1); 

Highfreqlim = High1 + High2;            
                                        
%% filter desing
% Normal frequencies and Amplitudes for filter desing
Niquist = find(fix(freq)<=24000,1,'last');
if ~isempty(High2);
    if Highfreqlim >= Niquist 
        ramp = linspace(PTF_i(High1+1),1,(Niquist-High1)/2);
        Amplitudes = [ones(Lowfreqlim-1,1);PTF_i(Lowfreqlim:High1);...
            ramp'; ones(Niquist-High1-length(ramp),1)];
    else
        Amplitudes = [ones(Lowfreqlim-1,1);PTF_i(Lowfreqlim:Highfreqlim);...
            ones(Niquist-Highfreqlim,1)];
    end
elseif isempty(High2);
    ramp = linspace(PTF_i(High1+1),1,(Niquist-High1)/2);
    Amplitudes = [ones(Lowfreqlim-1,1);PTF_i(Lowfreqlim:High1);...
        ramp'; ones(Niquist-High1-length(ramp),1)];
end
NormFreqs = freq(1:Niquist)./fs.*2;
% Arbitrary magnitude filter desing
h = fdesign.arbmag('N,F,A', FiltOrder, NormFreqs, Amplitudes);

% filter structure
Hd = design(h, 'freqsamp');

%% Plot time and frequency response of the resulting filters. 

InvFilt = freqz(Hd,nfft,'whole',fs);
%InvImpulseResp = Hd.Numerator;

fid = figure;
ax1 = axes('parent',fid,'position',[0.1 0.05 0.8 0.45]);
set(ax1,'box','on','linewidth',2, 'xgrid','on','ygrid','on',...
    'xlim',[20 fs/2],'ylim',[-30 30],'Fontsize',14);
ax=[20 fs/2 -30 30];
nr_of_dek=log10(ax(2)/ax(1));
dB_per_dek=(ax(4)-ax(3))/nr_of_dek;
whished_ratio=50;  %dB_per_dek.
xlabel(ax1,'Frequency [Hz]','Fontsize',14)
ylabel(ax1,'Gain [dB]','FontSize',14)
set(ax1,'PlotBoxAspectRatio',[whished_ratio/dB_per_dek 1 1],...
    'XTickLabel',{'100' '1000' '10000'},'XScale','log','Box','on');

l1 = line(freq,20*log10(abs(Ptf)),'color','g','LineWidth',2,'parent',ax1);
l2 = line(freq,20*log10(abs(Ptf_i)),'color','r','linewidth',2,'parent',ax1);
l3 = line(freq,20*log10(abs(InvFilt)),'color','b','linewidth',2,'parent',ax1);

legend(ax1,[l1 l2 l3],{'Measured Transfer function' 'Ideal Inverse'...
    'Designed filter'},'location','southwest')

tt = 1/fs:1/fs:length(Hd.Numerator)/fs;
ax2 = axes('parent',fid,'position',[0.1 0.6 0.8 0.32]);
set(ax2,'box','on','linewidth',2, 'xgrid','on','ygrid','on',...
    'xlim',[0 tt(end)],'ylim',[-0.7 1.5],'Fontsize',14);
xlabel(ax2,'Time [s]','Fontsize',14)
ylabel(ax2,'Amplitude','FontSize',14)

l4 = line(t,real(ptf),'color','g','LineWidth',2,'parent',ax2);
%l5 = line(t,ptf_m,'color','r','LineWidth',2,'parent',ax2);

l6 = line(tt,Hd.Numerator,'color','b','LineWidth',2,'parent',ax2);

legend(ax2,[l4 l6],{'Measured Impulse Response' ...
    'Inverse Filter'})
