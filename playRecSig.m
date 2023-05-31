%% function MicSig = playRecSig(auio,signal)
% Callback function for audioPlayerRecorder to play and record loger
% segments this is limited to variable sizes in memory. 
% TODO: read & write to files insted of memory.

function MicSig = playRecSig(auio,signal,CLS)
    % reshape the signal so that i fits into the buffer size of the
    % sound card
    if CLS
        pnoise = dsp.ColoredNoise(1,auio.BufferSize,1,'BoundedOutput',true);    
    end
    bf = auio.BufferSize;
    [sl, nch] = size(signal);
    FpB = floor(sl/bf);
    
    for i = 1:nch
        s1 = reshape(signal(1:FpB*bf,i),bf,FpB);
        % if there are samples that dont fill the last frame
        if rem(sl,bf) > 0
            Zpad = bf - rem(sl,bf);
            s2 = [signal(FpB*bf+1:end,i); zeros(Zpad,1)];
            AudioToDev(:,:,i) = [s1 s2];
        else 
            AudioToDev(:,:,i) = s1;
        end
    end
    % play each audio frames.
    InSignal = nan(size(AudioToDev));
    for i = 1:size(AudioToDev,2)
        if CLS
            NoiseSig = pnoise();
            x_ = auio([AudioToDev(:,i,1) 0.6*NoiseSig(:)]);
        else
            x_ = auio(squeeze(AudioToDev(:,i,:)));
        end
        for j = 1:nch
            InSignal(:,i,j) = x_(:,j);
        end
    end 
    MicSig = reshape(InSignal,size(InSignal,1)*size(InSignal,2),size(InSignal,3));
    
end