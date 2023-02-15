function cali_check_probeincalibrator(Ms,IS,fs)
    
                            display('------ CHECK PROBE SENSITIVITY IN PRESSURE CALIBRATOR ---------------');
                            display('Put ER-10C probe sound pressure calibrator B&K type 4130/4131.');
                            display('Press Enter.');
                            pause
    dur = 2;
    out = playrec('rec', 2*fs, 2);
    playrec('block',out);
    Out = playrec('getRec',out);

    rmsdig = sqrt(mean(Out(fs:end).^2));
    SPL = 20*log10(rmsdig/IS/Ms/20e-6);
    display(sprintf('Measured sound pressure is: %2.1f\n', roundn(SPL,-1)));
    
    %out = input('Type 0 if the levels are OK (default), 1 if not: ');
                            
    display('Sensitivity of probe microphone checked succesfully.');
    display('-------------------------------------------------------------');
end