function [repro RMS noiseRMS SNR] = teoae_band_processing(A,B,Fc)

%addpath /home/rop/Documents/Data/work/experiments/functions/octave/
%[b,a] = octdsgn(Fc,48e3,3);

Fs = 48000;
N = 16; % filter order
oneOctaveFilter = octaveFilter('FilterOrder', N, ...
    'CenterFrequency', Fc, 'Bandwidth', '1 octave',...
    'SampleRate', Fs);
OctFilt = getFilter(oneOctaveFilter);
A_av=OctFilt(mean(A,2));
B_av=OctFilt(mean(B,2));




repro=abs(corr(A_av,B_av))*100; % percentage

A_av=A_av'; B_av=B_av';
mean_response=(A_av+B_av)/2;
noise=(A_av-B_av)/2;

r=find(mean_response>0);
RMS=20*log10(sqrt(mean(mean_response(r(1):end).^2))/20e-6);
noiseRMS=20*log10(sqrt(mean(noise(r(1):end).^2))/20e-6);

SNR = RMS-noiseRMS;

end