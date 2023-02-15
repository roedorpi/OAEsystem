function dpdataprocess(subjectstring)

col = {'b-' 'r-' 'g-' 'm-'}; 
coll = {'bo-' 'rd-' 'gv-' 'm^-'}; 
folder = '../../DATA/';
dur = 1.3;
fs = 48e3;
num_reps = 5;
num_points = 16;

folder_contents = dir([folder '*.mat']);
%subjectstring = '^F04 DP R.*.mat$';

num_files = length(folder_contents);

k = 0;
Fig = figure;
hold on
axis([700 3500 -20 25])
box on
grid on
for i = 1:1:num_files
    filename = folder_contents(i).name;
    if regexp(filename, subjectstring) == 1
        ID = filename(1:3);
        k = k+1;
        load([folder filename]);
        DPLevel(:,:,k) = Data.DPLevel;
        DPFreq(:,k) = Data.DPFrequency;
        DPmean(:,k) = mean(DPLevel(:,:,k),1);
        DPstd(:,k) = std(DPLevel(:,:,k));
        area(DPFreq(:,k),Data.NoiseLevel,-20,'facecolor',[0.7 0.7 0.7])
        plot(DPFreq(:,k),DPmean(:,k), ...
            char(coll(k)),'linewidth',2,'markersize',8)
        errorbar(DPFreq(:,k),DPmean(:,k), ...
            DPstd(:,k),char(col(k)))
        title(['Subect: ' ID])
        
    end
end

[H,P,CI,STATS] = ttest(DPLevel(:,:,1),DPLevel(:,:,end),0.05)


% clearvars -except dur fs num_points num_reps timesignals
% 
% freqs = 0:1/dur:fs/2;
% freqdomain = 20*log10(2*abs(fft(timesignals'))/(dur*fs)/20e-6);
% freqdomain = freqdomain(1:dur*fs/2+1,:);
% size(freqdomain)
% semilogx(freqs,freqdomain(:,2))
% xlim([1 fs/2])
% grid on





