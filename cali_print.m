function cali_print()

    % LOADING
    load('./calibration.mat')
    Ms = mic_sens;
    IS = scIN;
    OSs = [scLOUT scROUT];
    LSs = [speakerL speakerR];
    
    display(['Ms  = ',num2str(Ms),' V/Pa']);
    display(['IS  = ',num2str(IS),' dig/V']);
    display(['OSL = ',num2str(OSs(1)),' V/dig']);
    display(['OSR = ',num2str(OSs(2)),' V/dig']);
    display(['LSL = ',num2str(LSs(1)),' Pa/V']);
    display(['LSR = ',num2str(LSs(2)),' Pa/V']);
    
end