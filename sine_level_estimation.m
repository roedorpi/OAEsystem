function res = sine_level_estimation(sig,fest,dur,fs)
    % sine_level_estimation uses goertzel to estimate 
    % the dB SPL level of one or more frequencies in 
    % a specified signal.

    % fs = 48e3;

    % f1 = (4250/1160).^((0:1:15)./(15))*1160;
    % f2 = f1*1.223;
    % fest = f1;%2006;
    % % t = 0:1/fs:3-1/fs;
    % % sig = sin(2*pi*fest*t);
    % sig = a;%e;

    SIG  = goertzel(double(sig),round(fest/fs*length(sig))+1);

    norm = dur*fs;%2*fs;%length(sig);
    res = 20*log10(2*abs(SIG)/norm/2e-5);
    % [roundn(f1',0) roundn(res,0)]


end