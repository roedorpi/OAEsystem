function calibration(~,~)
global calibrationFigure mainguiFigure calibration01 calibration02 calibration03

%% In/out System Calibration of computer vs. Edirol
%
%       SC = Edirol Sound Card
%       Do not use XLR cable as input
%
% (1)   Measurements of Probe Microphone sensitivity when the preamp is
%       adjusted to +20dB. Connect calibrator (1kHz 1Pa 94dB), and read the
%       measured RMS voltage on a multimeter.
%
%       If the measured value is 500 mVrms, everything is perfect. If it
%       differs, the next step should be scaled according to the error.
%
% (2)   Connect 1Vrms (If not an error has occured in (1)) 1kHz tone to the
%       input of the SC. 1Vrms is the maximum input which will be received 
%       for a preamp gain setting on +20dB, because only pressure levels up 
%       to 2Pa (100 dB) are measured undistorted. Adjust the input 
%       sensitivity to the point just before it peaks.
%
% (3)   Read out digital values on the computer and calculate the RMS
%       value (The first second of recorded data is skipped because of
%       peaks in the beginning).
%
% (4)   Disconnect the function generator, and measure the output value on 
%       the SC in Vrms, when a sinusoid of 1-kHz with an digital amplitude
%       of 1 is applied (Corrosponds to an digital rms amplitude of 
%       1/sqrt(2) = 0.707).
%
% (5)   Loop-back mode. Connect the input to the output of the SC. Run a
%       MLS sequence to get the impulse response and transfer function of 
%       the SC.
%
%% Nested functions
%
%   Playrec
%   AutoPlayrecInit
%   fastsine
%   fastmls
%   mlsgen
%

close(mainguiFigure);
scrsize=get(0,'ScreenSize');
calibrationFigure=figure;
set(calibrationFigure,'OuterPosition',[scrsize(3)/3  scrsize(4)/2.5 scrsize(3)/2.5 scrsize(4)/2],'color',[84 84 84]/256)
AutoPlayrecInit;

texte1 = { ['(1) Connect Edirol UA-25 Sound Card (SC) to computer']...
     ['(2) Turn the output to max on SC'] ['(3) If Mac, turn the volume on the computer to max']...
    ['(4) Push "MON SW" on SC'] ['(5) Adjust the sensitivity on SC']};

calibration01 = uicontrol(calibrationFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.8 0.8 0.15],...
    'BackgroundColor',[84 84 84]/256,...
    'ForegroundColor','w',...
    'FontSize',20,...
    'String','Before beginning:');
calibration02 = uicontrol(calibrationFigure,...
    'Style','text',...
    'Units','normalized',...
    'Position',[0.1 0.3 0.8 0.4],...
    'BackgroundColor',[84 84 84]/256,...
    'ForegroundColor','w',...
    'FontSize',12,...
    'String',texte1);

calibration03 = uicontrol(calibrationFigure,...
    'Style','pushbutton',...
    'Units','normalized',...
    'Position',[0.25 0.1 0.5 0.15],...
    'Fontsize',13,...
    'FontWeight','bold',...
    'String','Continue',...  
    'Callback',@cali_wrap);  
end