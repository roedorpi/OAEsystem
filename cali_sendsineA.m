function vals = cali_sendsineA(A,f,dur,fs)
%     dur = 10;
%     f = 500;
    

    t = 0:1/fs:dur-1/fs;
    x = A*sin(2*pi*f*t);
    
    input('Press Enter to send sines to both loudspeakers.')
    
    vals = zeros(1,2);
    for channel = 1:1:2
        pn = playrec('play',x',channel);
        playrec('block',pn);
        display(['Played 5 second in channel ',num2str(channel),'.']);
        vals(channel) = input('Vrms: ');
        if channel == 1
            display('Press Enter to send on channel 2');
            pause
        end
    end
    display('Sent sines to loudspeakers succesfully.')
    
end
