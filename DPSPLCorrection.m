function [dbspl] = DPSPLCorrection(DBSPL,FREQ,CHANNEL)

% Function that subtracts the transfer function of the probe speakers
% measured in the ear simulator. It is scaled so it has approximately 0 dB
% @1kHz. It can take 1 og 2 inputs for both DBSPL and FREQ. The mat-file 
% "DPSPLCorretionData.mat" must be included in the path of the main
% program. 
% CHANNEL = 1 => speaker 1 
% CHANNEL = 2 => speaker 2 

load DPSPLCorretionData;

if      CHANNEL == 1
        index1 = find(F >= FREQ(1),1);
        dbspl(1) = DBSPL(1) - H1(index1);
elseif  CHANNEL == 2
        index1 = find(F >= FREQ(1),1);
        dbspl(1) = DBSPL(1) - H2(index1);
end

M = size(DBSPL);
if      M(2) == 2 && CHANNEL == 1
        index2 = find(F >= FREQ(2),1);
        dbspl(2) = DBSPL(2) - H2(index2);
elseif  M(2) == 2 && CHANNEL == 2
        index2 = find(F >= FREQ(2),1);
        dbspl(2) = DBSPL(2) - H1(index2);
end

end
%% Used to generate data to this function
% clc; clear all; close all;
% 
% fs = 48000;
% 
% addpath /Users/martinkynde/Desktop/8_SEM_EE/8_Sem_SVN/Misc/Anders/lab/log
% 
% load('mls-1rep-4avg-16-480-coupler-LS1L-LS1LMIC-h_aver1isCMic-h_aver2isPMic.mat','h_aver1');
% ls1 = h_aver1;
% ls1 = 23.7411.*ls1; % Scaling factor applied to adjust ls1 to 0dB at 1 kHz
% 
% load('mls-1rep-4avg-16-480-coupler-LS2R-LS2RMIC-h_aver1isCMic-h_aver2isPMic.mat','h_aver1');
% ls2 = h_aver1;
% ls2 = 24.5.*ls2;    % Scaling factor applied to adjust ls2 to 0dB at 1 kHz
% 
% nfft = length(ls1);
% 
% % [xx,freq] = FREQZ(1,1,nfft,'whole',fs);
% 
% [H1,F] = freqz(ls1,1,nfft,fs);
% H1 = mag2db(abs(H1));
% 
% [H2,F] = freqz(ls2,1,nfft,fs);
% H2 = mag2db(abs(H2));
% 
% % Remove datapoints above 10 kHz
% index = find(F >= 10e3,1);
% F = F(1:index);
% H1 = H1(1:index);
% H2 = H2(1:index);
% 
% save('DPSPLCorretionData.mat','F','H1','H2');






