function Out = ProcessTEdata(H,File)

%%

%load([H.DataDir File]);


te.ID = File(1:3);
te.Ear = File(8);
Config = H.fio.Config;
Data = H.fio.Data;
MeasTime = datenum(regexp(File,...
        '\d\d-\d\d-\d\d\d\d\ \d\d\-\d\d-\d\d','match'),...
        'dd-mm-yyyy HH-MM-SS');
te.Mt = num2cell(datestr(MeasTime),2);
% Band-pass filter 
n = 6; Wn = [500 6000]*2/Config.Fs;
ftype = 'bandpass';
[b,a] = butter(n,Wn,ftype);

A = filtfilt(b,a,Data.A);
B = filtfilt(b,a,Data.B);

% 
A_av=mean(A,2);
B_av=mean(B,2);

[te.ABCorr te.RMS te.noiseRMS te.spectrum ...
    te.noise_spectrum te.mean te.noise te.Aav te.Bav] = teoae_processing(A,B);    
[te.AB_f_Corr(1) te.RMS_f(1) te.noiseRMS_f(1) te.SNR_f(1)] ...
    = teoae_band_processing(A,B,1000);
[te.AB_f_Corr(2) te.RMS_f(2) te.noiseRMS_f(2) te.SNR_f(2)] ...
    = teoae_band_processing(A,B,2000);
[te.AB_f_Corr(3) te.RMS_f(3) te.noiseRMS_f(3) te.SNR_f(3)] ...
    = teoae_band_processing(A,B,4000);


te.Fs = Config.Fs;
te.T = 1/te.Fs:1/te.Fs:length(te.mean)/te.Fs;
te.Freq = (0:length(te.mean)/2-1)/length(te.mean).*te.Fs;
te.A = A;
te.B = B;

Out.TEOAE_data = te;
end

