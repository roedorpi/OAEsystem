function fastsine(fs,channel)
% channel:    1 = left channel and 2 = right channel

%% Generating sinewave
amplitude = 1.0;            % Amplitude. Never above 1, else distortion!
T = 1/fs;                   % sampling period
t = 0:T:5-T;               % time vector
f = 1e3;

IN = amplitude*sin(2*pi*f*t);   % sinusoid

% Play
playBuffer = IN';

[pageNumber] = playrec('play',playBuffer, channel);

end