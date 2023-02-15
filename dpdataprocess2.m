function dpdataprocess(subjectstring)

col = {'b-' 'r-' 'g-' 'm-' 'c-' 'y-' 'b-' 'r-' 'g-' 'm-' 'c-' 'y-'}; 
coll = {'bo-' 'rd-' 'gv-' 'm^-' 'c<-' 'y>-' 'bo-' 'rd-' 'gv-' 'm^-' 'c<-' 'y>-'};
colll = {'b' 'r' 'g' 'm' 'c' 'y' 'b' 'r' 'g' 'm' 'c' 'y' }; 
folder = '../../DATA/';
dur = 1.3;
fs = 48e3;
num_reps = 5;

num_points = 16;

folder_contents = dir([folder '*.mat']);
num_files = length(folder_contents);
for i = 1:num_files
    aa(i) = folder_contents(i).datenum;
end
[aa indx] = sort(aa);
folder_contents = folder_contents(indx);

% f3 = figure(3);
% f3ax = axes('parent',f3);


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
        Ear{k} = filename(8);
        MeasTime(k) =  datenum(regexp(filename,...
        '\d\d-\d\d-\d\d\d\d\ \d\d\-\d\d-\d\d','match'),...
        'dd-mm-yyyy HH-MM-SS');
        load([folder filename]);
        [DPLevel(:,:,k) PriLev(:,:,k) FreQs(:,:,k)] = dpcalc([folder filename]);
        %DPLevel(:,:,k) = Data.DPLevel;
        DPFreq(:,k) = Data.DPFrequency;
        %DPmean(:,k) = mean(DPLevel(:,:,k),1);
        %DPstd(:,k) = std(DPLevel(:,:,k));
        Noise = cat(2, DPLevel(:,3,k)-DPLevel(:,4,k),DPLevel(:,3,k), DPLevel(:,3,k) + DPLevel(:,4,k));
        Y1 = cat(1,Noise(:,3),Noise(end:-1:1,2));
        Y2 = cat(1,Noise(:,1),Noise(end:-1:1,2));
        X = cat(1,DPFreq(:,k),DPFreq(end:-1:1,k));
        Z = ones(32,1);
        p1 = patch(X,Y1,Z);
        p2 = patch(X,Y2,Z);
        set(p1,'facecolor',char(colll(k)),'Edgecolor','none', 'facealpha',0.5)
        set(p2,'facecolor',char(colll(k)),'Edgecolor','none', 'facealpha',0.5)
        
        %plot(DPFreq(:,k),DPNoise(:,1,k),'k-','linewidth',1)
        plot(DPFreq(:,k),DPLevel(:,1,k), ...
            char(coll(k)),'linewidth',2, ...
            'markersize',8,'markerfacecolor', char(colll(k)))
        plot([2500 2600],[24-(k*1.5) 24-(k*1.5)],...
            char(col(k)),'linewidth',2)
        text(2610,24-(k*1.5), strcat(Ear(k),', ',datestr(MeasTime(k))))
        errorbar(DPFreq(:,k),DPLevel(:,1,k), ...
            DPLevel(:,2,k),char(col(k)))
        title(['Subect: ' ID])
    end
end

[H,P,CI,STATS] = ttest(DPLevel(:,:,1),DPLevel(:,:,end),0.05)

FileName = strcat('/media/data/oldnux_home/10.04_64/Desktop/SEC/figs/',... 
    handles.MainFigLegend);
PrintFig = figure;
ax = axes('parent',PrintFig);
ObjectHandles = get(handles.MainFig,'children');
PrintObjects = copyobj(ObjectHandles,ax);
set(ax,'xlim',[700 3500],'ylim',[-30 25], 'ytick', [-25:5:25], ...
    'box','on', 'Xgrid','on', 'ygrid','on', 'nextplot','add','fontsize',14);
%AX = [get(ax,'Xlim') get(ax,'ylim')];

xlabel('Frequency [Hz]','fontsize',14);
print(PrintFig,'-depsc',FileName);
close(PrintFig)
guidata(hObject, handles);
% clearvars -except dur fs num_points num_reps timesignals
% 
% freqs = 0:1/dur:fs/2;
% freqdomain = 20*log10(2*abs(fft(timesignals'))/(dur*fs)/20e-6);
% freqdomain = freqdomain(1:dur*fs/2+1,:);
% size(freqdomain)
% semilogx(freqs,freqdomain(:,2))
% xlim([1 fs/2])
% grid on





