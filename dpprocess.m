function dpprocess(filename, handles)

col = {'b-' 'r-' 'g-' 'm-' 'c-' 'y-'}; 
coll = {'bo-' 'rd-' 'gv-' 'm^-' 'c<-' 'y>-'};
colll = {'b' 'r' 'g' 'm' 'c' 'y'}; 

k = handles.numberofplots;
axes(handles.MainFig)
hold on
axis([700 3500 -20 25])
box on
grid on
ID = filename(1:3);
Ear = filename(8);
MeasTime =  datenum(regexp(filename,...
    '\d\d-\d\d-\d\d\d\d\ \d\d\-\d\d-\d\d','match'),'dd-mm-yyyy HH-MM-SS');

[DPLevel PriLev FreQs TS TS_f fs T] = dpcalc(filename);
DPFreq = FreQs(:,3);
Noise = cat(2, DPLevel(:,3)-DPLevel(:,4),DPLevel(:,3), DPLevel(:,3) + DPLevel(:,4));
Y1 = cat(1,Noise(:,3),Noise(end:-1:1,2));
Y2 = cat(1,Noise(:,1),Noise(end:-1:1,2));
X = cat(1,DPFreq,DPFreq(end:-1:1));
Z = ones(2*length(DPFreq),1);
p1 = patch(X,Y1,Z);
p2 = patch(X,Y2,Z);
set(p1,'facecolor',char(colll(k)),'Edgecolor','none', 'facealpha',0.5)
set(p2,'facecolor',char(colll(k)),'Edgecolor','none', 'facealpha',0.5)

plot(DPFreq,DPLevel(:,1),char(coll(k)),'linewidth',2, ...
            'markersize',8,'markerfacecolor', char(colll(k)))
plot([2500 2600],[24-(k*1.5) 24-(k*1.5)],char(col(k)),'linewidth',2)
text(2610,24-(k*1.5), strcat(ID, Ear,', ',datestr(MeasTime)))
errorbar(DPFreq,DPLevel(:,1), DPLevel(:,2),char(col(k)))
end


%[H,P,CI,STATS] = ttest(DPLevel(:,:,1),DPLevel(:,:,end),0.05)


% clearvars -except dur fs num_points num_reps timesignals
% 
% freqs = 0:1/dur:fs/2;
% freqdomain = 20*log10(2*abs(fft(timesignals'))/(dur*fs)/20e-6);
% freqdomain = freqdomain(1:dur*fs/2+1,:);
% size(freqdomain)
% semilogx(freqs,freqdomain(:,2))
% xlim([1 fs/2])
% grid on





