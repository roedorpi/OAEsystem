function IS = cali_check_input(IS,Ms,fs)
         
    
                            % Don't-care stuff indented
                            display('------------------------ CHECK INPUT ------------------------');
                            display('Put ER-10C probe in 1 kHz, 1 Pa (94 dB SPL) calibrator.');
                            display(['Using Ms = ',num2str(Ms),' V/Pa.']);
                            display('Press Enter.');
                            pause
                            display('Measuring sound pressure level...');

    pn = playrec('rec',2*fs,2);
    playrec('block',pn);
    x = playrec('getRec',pn);

    rmsdig = sqrt(mean(x(fs:end).^2));
    SPL = 20*log10(rmsdig/IS/Ms/20e-6);

                            display(['RESULTS.']);
                            display(['Measured:    ',num2str(SPL),' dB SPL.']);
                            display(['Measured IS: ',num2str(rmsdig/Ms),' dig/V.']);
                            choice = input('Type 0 to accept SPL (default), 1 to update IS: ')
                            if choice == 1
                                IS = rmsdig/Ms; 
                                display('Updated IS.')
                            else
                                display('Accepted SPL. IS not updated.')
                            end
                            display('Checked input succesfully.')
                            display('-------------------------------------------------------------');

end