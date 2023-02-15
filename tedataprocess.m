function tedataprocess(subjectstring)


folder = '../../DATA/';
dur = 1.3;
fs = 48e3;
num_reps = 512;
TimeRejection = 0.0038;

folder_contents = dir([folder '*.mat']);
%subjectstring = '^F08 TE R.*.mat$';

num_files = length(folder_contents);

k = 0;
Fig = figure;
hold on
xlabel('Time [s]');
ylabel('Emission Amplitud [Pa]');
%axis([700 3500 -20 25])
box on
grid on
% dh = fdesign.highpass(50*2/fs,500*2/fs,50,0.001);
% dl = fdesign.lowpass(6e3*2/fs,8e3*2/fs,0.001,60);
% Hh = design(dh,'equiripple');
% Hl = design(dl,'equiripple');
% hgd = grpdelay(Hh); hgd = round(hgd(1));
% lgd = grpdelay(Hl); lgd = round(lgd(1));
% hlgd = hgd + lgd;
% Himp = impz(Hh);
% Limp = impz(Hl);

n = 6; Wn = [500 6000]*2/fs;
ftype = 'bandpass';
[b,a] = butter(n,Wn,ftype);
h1=dfilt.df2(b,a);   

%[Num Den] = butter(10,12e3/fs);
for i = 1:1:num_files
    filename = folder_contents(i).name;
    if regexp(filename, subjectstring) == 1
        ID = filename(1:3);
        k = k+1;
        load([folder filename]);
        A(:,k) = mean(filtfilt(b,a,Data.A),2);
        B(:,k) = mean(filtfilt(b,a,Data.B),2);
        Corr = corrcoef(A(:,k),B(:,k));
        ABCorr(:,k) = Corr(1,2);
        AB(:,k) = mean([A(:,k), B(:,k)],2);
        t = 1/fs:1/fs:size(A,1)/fs;
        leg{k} = sprintf('AB correrlation: %1.4f', ABCorr(k));
    end
end
%plot(t(1:end-hlgd),AB(hlgd+1:end,:))
plot(t,AB)
legend(leg);
title(['Subect: ' ID])
%Corr_reppetitions = corrcoef(AB)





