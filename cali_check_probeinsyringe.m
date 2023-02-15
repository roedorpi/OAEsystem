function cali_check_probeinsyringe(OSs,Ms,IS,reference,fs)
    
                            display('--------------- CHECK PROBE LEVELS IN SYRINGE ---------------');
                            display('Put ER-10C probe in 1cm^3 syringe.');
                            display('Press Enter.');
                            pause
    dur = 2;
    f = 500;
    AVrms = 1;

    A = AVrms*sqrt(2);
    t = 0:1/fs:dur-1/fs;
    for channel = 1:1:2
        x = A/OSs(channel)*sin(2*pi*f*t);
                            display(['Playing and recording a ',num2str(AVrms),'Vrms@',num2str(f),'Hz tone from channel ',num2str(channel),'.'])
        pn = playrec('playrec',x',channel,dur*fs,2);
        playrec('block',pn);
        y = playrec('getRec',pn);
        y = y(fs:end)/IS/Ms;
        rmsydB = 20*log10(sqrt(mean(y.^2))/20e-6);
                            display(['Reference: ',num2str(reference(channel)),...
                                     ' dBSPL, Measured: ',num2str(rmsydB),...
                                     ', Difference (R-M): ',num2str(reference(channel)-rmsydB),'.']);
    end
                       
                            display('Checked probe levels in syringe succesfully.');
                            display('-------------------------------------------------------------');

end
