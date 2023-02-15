function vals = cali_sendsineVrms(Vrms,f,dur,OSs,fs)

    
    input('Press Enter to send sines to both loudspeakers.')
    

    t = 0:1/fs:dur-1/fs;

    vals = zeros(1,2);
    for channel = 1:1:2
        x = Vrms*sqrt(2)/OSs(channel)*sin(2*pi*f*t);
        pn = playrec('play',x',channel);
        playrec('block',pn);
        display(['Played 5 second in channel ',num2str(channel),'.']);
        vals(channel) = input('dB SPL: ');
        if channel == 1
            display('Press Enter to send on channel 2');
            pause
        end
    end
    display('Sent sines to loudspeakers succesfully.')
    
end
