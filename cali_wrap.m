function cali_wrap(varargin)%doSy,doMs,doIS,doOS,doLS)
    
    doSy = 0; doMs = 0; doIS = 0; doOS = 0; doLS = 0; doCK = 0;

    args = varargin;
    for i = 1:1:length(args)
        doSy = xor(doSy,strcmp(args(i), 'syringe'));
        doMs = xor(doMs,strcmp(args(i), 'micsens'));
        doIS = xor(doIS,strcmp(args(i), 'inputscaling'));
        doOS = xor(doOS,strcmp(args(i), 'outputscaling'));
        doLS = xor(doLS,strcmp(args(i), 'lssens'));
        doCK = xor(doCK,strcmp(args(i), 'miccheck'));
    end
%     [doSy doMs doIS doOS doLS] % For testing


    fs = 48e3;                  % Temporary hard-coding
    AutoPlayrecInit(fs);
    
%     AutoPlayrecInit;

    % LOADING
    load('./calibration.mat')
    Ms = mic_sens;
    IS = scIN;
    OSs = [scLOUT scROUT];
    LSs = [speakerL speakerR];
    
    refinSyringe = [78.9169 78.3397]; % dB SPL measured in syringe, when 
                                  %   levels are correct in ear simulator.

    %%%%%%%%%
    display('Loaded conversion factors from calibration file.')
    
    cali_print();
    
    
    if doSy == 1
        
        cali_check_probeinsyringe(OSs,Ms,IS,refinSyringe,fs);
        display('Info: This code does not react on the response.')
        
    end
    
    if doMs == 1
        display('MEASURE Ms')
        display('Put probe in 1Pa@1kHz calibrator and connect preamp to voltmeter.')
        choice = input('Measured RMS voltage (-1 to cancel): ');
        if choice == -1
            display('Not using newly measured Ms.');
        else
            Ms = choice;
            display('Using newly measured Ms.')
        end
    end
    
    if doIS == 1
        IS = cali_check_input(IS,Ms,fs);
    end
    
    if doOS == 1
        % OSs:
        % Connect output to volt meter.
        display('Measure new OSs:');
        display('Connect output to volt meter.');
        OSsnew = cali_sendsineA(1,500,5,fs)*sqrt(2);
        display('Changed OSs.');
        
        choice = input('Type 0 to keep old OSs (default), 1 to save.');

        if choice == 1
            OSsold = OSs;
            OSs = OSsnew;
        end
    end
    
    if doLS == 1
        % LSs:
        % Put probe in ear simulator attached to SPL meter.
        display('Measure new LSs:');
        display('Put probe in ear simulator attached to SPL meter.');
        if exist('OSsold','var') == 0
            cOSs = OSs;
        else
            choice = input('Type 0 (default) to use new OSs just measured, 1 to use old: ');
            if choice == 1
                cOSs = OSsold;
            else
                cOSs = OSs;
            end
        end
        LSsnew = 20e-6*10.^(cali_sendsineVrms(1,500,5,cOSs,fs)/20);
        display('Changed LSs.');

        choice = input('Type 0 to keep old LSs (default), 1 to save.');

        if choice == 1
            LSs = LSsnew;
        end
    end
    
    
    if doCK == 1
        % Check system:
        % Put probe in pressure calibrator @1kHz 94 dB SPL;.
        cali_check_probeincalibrator(Ms,IS,fs);
        display('Info: This code does not react on the response.')
        
    end
    
    playrec('reset')
    
    
    % SAVING
%     mic_sens = Ms;
%     scIN = IS;
%     scLOUT = OSs(1);
%     scROUT = OSs(2);
%     speakerL = LSs(1);
%     speakerR = LSs(2);
% 
%     save('./calibration','mic_sens','scIN',...
%                                      'scLOUT','scROUT',...
%                                      'speakerL','speakerR');
%     display('Saved new conversion factors as specified.')
%     
%     cali_print();
%     %%%%%%%%
end