function y = bandpass_denoise(x,Fpass1,Fpass2, Astop1, Astop2, Fs)
%UNTITLED Returns a discrete-time filter object.

%
% MATLAB Code
% Generated by MATLAB(R) 7.11 and the Signal Processing Toolbox 6.14.
%
% Generated on: 15-Aug-2011 15:04:39
%

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
%Fpass1 = 500;                            % First Passband Frequency
Fstop1 = Fpass1 - Fpass1^(2/3);         % First Stopband Frequency
%Fpass2 = 9000;        % Second Passband Frequency
Fstop2 = Fpass2 + Fpass2^(2/3);       % Second Stopband Frequency
%Astop1 = 60;          % First Stopband Attenuation (dB)
Apass  = 1;           % Passband Ripple (dB)
%Astop2 = 80;          % Second Stopband Attenuation (dB)
match  = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2, Fs);
Hd = design(h, 'butter', 'MatchExactly', match);

y = filter(Hd,x);

% [EOF]
