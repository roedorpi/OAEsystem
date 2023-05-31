%%
impulses = [impulses, zeros(1,0.01*fs)];
%%
ax = axes(Parent=figure(2));
%%
delete(ax.Children);    
pl(1) = line('parent',ax,'xdata',1:length(impulses),'ydata',impulses,'color','b');
pl(2) = line('parent',ax,'xdata',1:length(MicIN(363:end,1)),'ydata',MicIN(363:end,1));
%%
TStart = tic;
T = 0;
while T < 10
    T = toc(TStart);
    MicIN = playRecSig(h.AuIO,[impulses',zeros(length(impulses),1)]);
    
    set(pl(2),"YData",-MicIN(363:end,1))

    drawnow
end