function Out = ProcessTE_CLSdata(fio)

%%
Measurements = who(fio);
o = 1;
p = 1;
for i = 1:length(Measurements)
    m = fio.(Measurements{i});
    Config = m.Config;
    Data = m.Data;
    if strcmp(Data.Ear,'L')
        Ear = 1;
    else
        Ear = 2;
    end
    te.Mt = [Data.Date ' ' Data.EndTime];
    % Band-pass filter 
    n = 6; Wn = [500 6000]*2/Config.Fs;
    ftype = 'bandpass';
    [b,a] = butter(n,Wn,ftype);
    Stim = filtfilt(b,a,Data.Stimulus);
    
    for j = 1:size(Stim,2)  
        if abs(max(Stim(:,j))) > abs(min(Stim(:,j)))
            Stim(:,j) = -Stim(:,j);
        end
    end
    Stim = reshape(Stim,20e-3*Config.Fs,4,[]);
    Stim(1:round(4.5e-3*Config.Fs),:,:) = 0;
    A = squeeze(mean(Stim(:,:,1:2:end),2));
    B = squeeze(mean(Stim(:,:,2:2:end),2));
    [te.ABCorr, te.RMS, te.noiseRMS, te.spectrum, ...
        te.noise_spectrum, te.mean, te.noise, te.Aav, te.Bav] = teoae_processing(A,B);    
    [te.AB_f_Corr(1), te.RMS_f(1), te.noiseRMS_f(1), te.SNR_f(1)] ...
        = teoae_band_processing(A,B,1000);
    [te.AB_f_Corr(2), te.RMS_f(2), te.noiseRMS_f(2), te.SNR_f(2)] ...
        = teoae_band_processing(A,B,2000);
    [te.AB_f_Corr(3), te.RMS_f(3), te.noiseRMS_f(3), te.SNR_f(3)] ...
        = teoae_band_processing(A,B,4000);
    te.T = 1/Config.Fs:1/Config.Fs:length(te.mean)/Config.Fs;
    te.Freq = (0:length(te.mean)/2-1)/length(te.mean).*Config.Fs;
    te.A = A;
    te.B = B;
    
    if Config.CLS
        Out.TEOAE_CLS(p,Ear) = te;
        p = p+1;
    else
        Out.TEOAE(o,Ear) = te;
        o = o+1;
    end
end

