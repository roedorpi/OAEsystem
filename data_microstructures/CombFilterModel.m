function varargout = CombFilterModel(varargin)

%%
% Analysis frequencies

D = varargin{1}; % subject data 
fa = varargin{2};  %[1 1.5 2 3 4 5]*10^3 %audiometric frequencies

%fc = varargin{3}; %[1001 1501 2002 3003 4004 5005]; % measured f2 frequencies
for i = 1:length(fa), fc(i) = D.f2(find(D.f2>=fa(i),1,'first')); end


% bandwidth for analysis, 2/32 (1),4/42 (2),4/32 (3),6/32 (4)')
for i = 1:length(varargin{3})
   switch varargin{3}(i)
       case 1
            BandWidth(i) = 2/32;
       case 2
           BandWidth(i) = 4/42;
       case 3
           BandWidth(i) = 4/32;
       case 4 
           BandWidth(i) = 6/32;
   end
end

G = 10^(3/10);

% Number of observations
NrObs = varargin{4};

% Phase estimate multiplier
PhMult = varargin{5};

%non-linear fit  and options

% ft = fittype(...
%     'sqrt((a+b)^2*cos(2*pi*x/343 - c/2).^2 + (a-b)^2*sin(2*pi*x/343 -c/2).^2);',...
%     'independent', 'x', 'dependent', 'y' );
% ft = fittype(...
%     'abs((a+b*cos(2*pi*x*c))-1i*(b*sin(2*pi*x*c)));',...
% %     'independent', 'x', 'dependent', 'y' );
% ft = fittype(...
%     'sqrt(a^2 + 2*a*b*cos(2*pi*x*c) + b^2);',...
%     'independent', 'x', 'dependent', 'y' );
ft = fittype(...
    'sqrt(a^2 + 2*a*b*cos(2*pi*x*c) + b^2);',...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions(ft); 

%opts.Method =  'NonlinearLeastSquares';
%opts.Algorithm = 'Levenberg-Marquardt';
%opts.Display = 'iter';
%opts.Robust = 'Bisquare';
opts.TolFun = 10e-16;
opts.TolX = 10e-16;
%% 
for i = 1:length(fa) % center frequencies for analysis 
    f1 = fc.*G^(-1/2*BandWidth(i)); %lower band edge
    f2 = fc.*G^(+1/2*BandWidth(i)); %upper band edge
    DeltaF = f2-f1;
    
    opts.Lower = [0 0 -pi];
    opts.Upper = [Inf Inf pi];
    f = floor(f1(i)):floor(f2(i));
   %fprintf('Freq: %i [Hz]\n',fa(i))
    % desired frequency points
    DesF = round(linspace(D.f2(find(D.f2>=f(1),1,'first')),...
        D.f2(find(D.f2<=f(end),1,'last')),NrObs(i)));
    MM = D.f2(find(D.f2>=f(1),1,'first'):find(D.f2<=f(end),1,'last'))';

    if length(MM) < NrObs(i)
    %    fprintf('Insufficient data\n')
    elseif length(MM) == NrObs(i)
        MF = find(D.f2>=f(1),1,'first'):find(D.f2<=f(end),1,'last');
    elseif length(MM) > NrObs(i)
         for l = 1:NrObs(i)
             FreqDiff = MM - DesF(l);
             if find(FreqDiff == 0)
                MF(l) = find(D.f2==MM(find(FreqDiff == 0,1,'first')),1,...
                    'first');
             else 
                [~,b] = min(abs(FreqDiff)); 
                MF(l) = find(D.f2==MM(b),1);
             end
         end
    end  
    if exist('MF','var')
        freq = D.f2(MF);
        DPdB = D.dpoae(MF); % dpoae in dB
        dp = 10.^(DPdB/20)*20e-6; % values in Pa.        
        opts.StartPoint = [mean(dp),0.5*mean(dp), 1/(DeltaF(i)*PhMult(i))];
        [Fit,Gof,Out] = fit(freq,dp,ft,opts);
        r2 = Gof.rsquare;
        Step = 0.01;
        PhaAdj = PhMult(i);
%         ModdB = 20*log10(feval(Fit,freq)./20e-6);
%         for j = 1:length(DPdB)-1
%             if abs(DPdB(j)-DPdB(j+1)) <= 1
%                 SlopeDP(j) = 1;
%             elseif DPdB(j) > DPdB(j+1)
%                 SlopeDP(j) = -1;
%             elseif DPdB(j) < DPdB(j+1) 
%                 SlopeDP(j) = 1;
%             end
%             if abs(ModdB(j)-ModdB(j+1)) <= 1
%                 SlopeMod(j) = 1;
%             elseif ModdB(j) > ModdB(j+1) 
%                 SlopeMod(j) = -1;
%             elseif ModdB(j) < ModdB(j+1) 
%                 SlopeMod(j) = 1;
%             end
%         end
%         %NumZeros = find(SlopeDP == 0);
        %NumOther = find(SlopeDP ~= 0);
        while r2 <= 0.9 &&  PhaAdj <= 3.1 % && PhaAdj >= 0.3 % && Resid > MaxErr
%             %if isempty(NumZeros)  % all differences greater that 1 dB
%                 if SlopeDP(1) == SlopeDP(3) && SlopeMod(1) ~= SlopeMod(3)
                    PhaAdj = PhaAdj + Step;
%                 elseif SlopeDP(1) ~= SlopeDP(3) && SlopeMod(1) == SlopeMod(3)   
%                     PhaAdj = PhaAdj + Step;
%                 else
%                     PhaAdj = PhaAdj/0.5;
%                 end
            %elseif numel(NumZeros) == 1 
            %    if NumZeros == 1 && sum(SlopeDP(NumOther))~= 0
                    %% no peak, nonotonic up or down    
            %    else  
                    %% one peak
            %    end
            %else
                %% no change 
            %end
            opts.StartPoint = [mean(dp),0.5*mean(dp), 1/(DeltaF(i)*PhaAdj)];
            [Fit,Gof,Out] = fit(freq,dp,ft,opts);
%             ModdB = 20*log10(feval(Fit,freq)./20e-6);
%             for j = 1:length(DPdB)-1
%                 if abs(ModdB(j)-ModdB(j+1)) <= 1
%                     SlopeMod(j) = 1;
%                 elseif ModdB(j) > ModdB(j+1)
%                     SlopeMod(j) = -1;
%                 elseif ModdB(j) < ModdB(j+1)
%                     SlopeMod(j) = 1;
%                 end
%             end
            r2 = Gof.rsquare;
            H = getappdata(gcf,'H');
            PhMult = str2num(H.ModParam(7).String);
            PhMult(i) = PhaAdj;
            H.ModParam(7).String = strcat('[',num2str(PhMult),']');
            r_2 = num2str(Gof.rsquare);
            rmse_ = num2str(Gof.rmse);
            H.Table.Data(4*(i-1)+1,3:4) = {r_2, rmse_}; 
            setappdata(H.mainfig,'H',H);
            
        end

        F = freq(1):freq(end);
        Mod.ModelEstimate{i} = 20*log10(Fit(F)/20e-6);
        Mod.ModelFreq{i} = F;
        Mod.DPEnergyAverage(i) = 20*log10(mean(dp)/20e-6); 
        Mod.DPMeasuredFrequency(i) = ...
            20*log10(dp(find(freq>=fa(i),1,'first'))/20e-6);
        Mod.DPFreq(i) = freq(find(freq>=fa(i),1,'first'));
        
        Mod.DPEstimGenCompAmp(i) = 20*log10(max([Fit.a Fit.b])/20e-6); %,subj) = 
        Mod.InputIndx{i} = MF;
        Mod.Fit{i} = Fit;
        Mod.Gof{i} = Gof;
        Mod.Out{i} = Out;

    else
        Mod.InputIndx{i} = [];
    end
    clear MF
        
end
varargout{1} = Mod;

