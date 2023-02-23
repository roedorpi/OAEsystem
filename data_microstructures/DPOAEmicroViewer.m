function varargout = DPOAEmicroViewer(varargin)
%% DPOAEmicroViewer
% GUI to inspect and fit a Comb filter model (CombFilter.m) to DPOAE micro 
% structures. All micro DPOAE data measured by Karen Reuter are collected 
% the data set loaded by the GUI.
% 
% How to use the GUI:
% In MATLAB go to the folder where the folowing files are saved: 
% DPmicro.mat : MATLAB data file, measured data from 93 subjects.
% CombFilter.m : Model definition
% CombFilterModel.m : Fitting of comb filter to the selected data.
% DPOAEmicroViewer.m : This file.
% Execute DPOAEmicroViewer.m either from the editor or the command window.
%
% Once the GUI is loaded, load data from one subject using the top dropdown
% menu to the left of the axes. This will load the chosen subject's data 
% and display the DPOAE measured with high frequency resolution, showing 
% both DPOAE level and estimated noise.
%
% Model Parameters:
% The panel named "Model Parameters" can be used to determine the amount of
% data and at which frequencies the analysis will be carried out. The
% following is a description of the parameters that are applied to the
% model when ever the file CombFilterModel.m is executed by pressing the
% "Calculate Model" button at the bottom of this panel.
% 
%  # Audiometric Frequencies: Standard audiometric frequencies to carry out 
%  the analysis
%  # 

load('DPmicro.mat');
H.D = D;
clear D


ScrSize = get(0,'ScreenSize');
H.mainfig = figure(...
    'Outerposition', ScrSize,...
    'name', 'DPOAE Micro Viewer',...  % change this to fit your experiment
    'resize', 'on',...
    'toolbar', 'auto',...
    'visible', 'on',...
    'tag', 'mainfig',...
    'MenuBar', 'figure',...
    'DockControls','on'); %     'MenuBar', 'none',...

    H.AX = axes('nextplot','add','ylim',[-20 20],...
        'xlim',[650 7000],'xscale','log','box','on',...
        'xtick',[1000 1500 2000 3000 4000 5000],'Xminortick','off',...
        'units','normalized','Position',[0.35,0.1 0.65 0.8]);
    ylabel(H.AX,'DPOAE amplitude [dB SPL]');
    xlabel(H.AX,'f_2 [Hz]');
    TITLE = sprintf('Subject: %i',[]);
    dBPerDec = 50;
    pbaspect(H.AX,[diff(H.AX.XLim)/(H.AX.XLim(1)*10), diff(H.AX.YLim)/dBPerDec 1])
    title(H.AX,TITLE);
    H.Subject = uicontrol('Parent',H.mainfig,'style','popup','units','normalized',...
        'Position', [0.075,0.7,0.175 0.2],'string',[{'Choose a Subject'},num2cell(1:93)],'callback',@Loadsubj);
        
    
    %% Model parameters
    H.ModEval= uipanel('Parent',H.mainfig,'units','normalized',...
        'Position', [0.025,0.62,0.3 0.25],'title','Model Parameters');
    % Frequencies
    H.ModTxt(1) = uicontrol('Parent',H.ModEval,'style','text','units','normalized',...
        'Position', [0.05,0.92,0.9 0.05],'string','Audiometric Frequencies [Hz]');
    H.ModParam(1) = uicontrol('Parent',H.ModEval,'style','edit','units','normalized',...
        'Position', [0.05,0.82,0.9 0.1],'string','[1000 1500 2000 3000 4000 5000]');   
    % bandwidth
    H.ModTxt(3) = uicontrol('Parent',H.ModEval,'style','text','units','normalized',...
        'Position', [0.05,0.66,0.4 0.15],'string','Bandwidth: 2/32 (1),4/42 (2),4/32 (3),6/32 (4)');
    H.ModParam(3) = uicontrol('Parent',H.ModEval,'style','edit','units','normalized',...
        'Position', [0.05,0.565,0.4 0.1],'string','[2 1 1 1 1 1]');
    % Number of observations
    H.ModTxt(4) = uicontrol('Parent',H.ModEval,'style','text','units','normalized',...
        'Position', [0.55,0.65,0.4 0.1],'string','Number of observations');
    H.ModParam(4) = uicontrol('Parent',H.ModEval,'style','edit','units','normalized',...
        'Position', [0.55,0.565,0.4 0.1],'string','[4 4 4 4 4 4]');
    % Phase multiplier 
    H.ModTxt(7) = uicontrol('Parent',H.ModEval,'style','text','units','normalized',...
        'Position', [0.05,0.27,0.9 0.2],'string',...
        'Initial phase estimate, taken as the bandwidth under analysis multiplied by the constant. One value for each frequency.');
    H.ModParam(7) = uicontrol('Parent',H.ModEval,'style','edit','units','normalized',...
        'Position', [0.05,0.18,0.9 0.1],'string','[0.5 0.5 0.5 0.5 0.5 0.5]');
    % Execute button
    H.ModParam(8) = uicontrol('Parent',H.ModEval,'style','pushbutton',...
         'units','normalized','Position', [0.05,0.01,0.9 0.15],...
         'string','Calculate Model','callback',@ModelCalculation,'enable','off');

     
 %% Regresion results
 %
    H.ModRes = uipanel('Parent',H.mainfig,'units','normalized',...
            'Position', [0.025,0.1,0.3 0.525],'title','Regresion result');
    Freqs = [str2num(H.ModParam(1).String)]' ;
    m = 1;
    for i = 1:4:length(Freqs)*4-3
        TabFreqs{i} = num2str(Freqs(m));
        m = m + 1;
        TabFreqs{i + 1} = 'A [Pa]';
        TabFreqs{i + 2} = 'B [Pa]';
        TabFreqs{i + 3} = 'Ph. [Hz]';
        [TabData{i,1:2}] = deal('%%%%%');
        [TabData{i + 1,1:2}] = deal('--');
        [TabData{i + 2,1:2}] = deal('--');
        [TabData{i + 3,1:2}] = deal('--');
        [TabData{i,3:4}] = deal('--');
        [TabData{i + 1,3:4}] = deal('%%%%%');
        [TabData{i + 2,3:4}] = deal('%%%%%');
        [TabData{i + 3,3:4}] = deal('%%%%%');
    end

    H.Table = uitable(H.ModRes,'Data', TabData);     
    H.Table.ColumnFormat = {'char','char'};
    H.Table.Units = 'normalized';
    H.Table.Position = [0.01 0.01 0.98 0.98];

    H.Table.RowName = TabFreqs;%;
    H.Table.ColumnName = {'Coeff.','95% C.I. for Coeff.','r^2','RMSE'}; 
    H.Table.ColumnWidth = {90 140 80 80};     
    setappdata(H.mainfig,'H',H)
    varargout{1} = H;
end


function Loadsubj(hObject,~,~)
    H = getappdata(hObject.Parent,'H');
    delete([H.AX.Children]);
    H.Subj = hObject.Value -1;
    if H.Subj > 0 && H.Subj <= length(H.D)
        H.L(2) = line(H.D(H.Subj).f2,H.D(H.Subj).noise,'parent',H.AX,'color',[0.8 0.8 0.8],'linewidth', 2);      
        H.L(1) = line(H.D(H.Subj).f2,H.D(H.Subj).dpoae,'parent',H.AX,'color',[0.6 0.6 0.6],'linewidth', 2);
        H.AX.Title.String = sprintf('Subject: %i',H.Subj);
        H.Legend = legend(H.AX,H.L,{'dpoae','noise'},'Location','northwest');
        H.ModParam(8).Enable = 'on';
        H.ModParam(7).String = '[0.5 0.5 0.5 0.5 0.5 0.5]';
    end
    %set(H.Estimates.Children,'Enable','on')
    setappdata(H.mainfig,'H',H)
end

function ModelCalculation(hObject,~,~)
    H = getappdata(get(hObject.Parent,'Parent'),'H');
    % clear previous model
    if isfield(H, 'M')
        delete(H.M);
    end
    % clear previous model plot
    if isfield(H,'EstMarker')
        delete(H.EstMarker)
        delete(H.M)
        H.Legend = legend(H.AX,H.L,{'dpoae','noise'},'Location','northwest'); 
    end
    if H.Subj > 0
        H.fa = str2num(H.ModParam(1).String);
        BW = str2num(H.ModParam(3).String);
        NrObs = str2num(H.ModParam(4).String);
        CutoffFreq = str2num(H.ModParam(7).String);
        H.Mod = CombFilterModel(H.D(H.Subj),H.fa,BW,NrObs,CutoffFreq);
        
        
        for j = 1:length(H.fa)
            if ~isempty(H.Mod.InputIndx{j})
            H.M(j,2) = line(H.AX,H.D(H.Subj).f2(H.Mod.InputIndx{j}),...
                H.D(H.Subj).dpoae(H.Mod.InputIndx{j}),'markersize',6,...
            'Marker','o','Color',[0.1 0.6 0.3],...
            'markerfacecolor',[0.1 0.6 0.3],'linewidth', 2);
            H.M(j,1) = line(H.Mod.ModelFreq{j},...
                H.Mod.ModelEstimate{j},'parent',H.AX,...
                'color',[0.3 0.1 0.6],'linewidth', 2);
           
            %A = num2str(roundn(20*log10(abs(H.Mod.Fit{j}.a)/20e-6),-1));
            %B = num2str(roundn(20*log10(abs(H.Mod.Fit{j}.b)/20e-6),-1));
            % Coeficients
            A = num2str(roundn(H.Mod.Fit{j}.a,-6));
            B = num2str(roundn(H.Mod.Fit{j}.b,-6));
            Phase = num2str(roundn(1/H.Mod.Fit{j}.c,-1));
            H.Table.Data(4*(j-1)+2:4*(j-1)+4,1) = {A, B, Phase}'; 
            
            % confidence interval for coeficients
            CI = confint(H.Mod.Fit{j});
            A_ci = ['(',num2str(roundn(CI(1,1),-6)),';',num2str(roundn(CI(2,1),-6)),')'];
            B_ci = ['(',num2str(roundn(CI(1,2),-6)),';',num2str(roundn(CI(2,2),-6)),')'];
            Phase_ci = ['(',num2str(roundn(1/CI(1,3),-1)),';',num2str(roundn(1/CI(2,3),-1)),')'];
            H.Table.Data(4*(j-1)+2:4*(j-1)+4,2) = {A_ci, B_ci, Phase_ci}'; 
            
            % Goodness of fit
            r2 = num2str(H.Mod.Gof{j}.rsquare);
            rmse = num2str(H.Mod.Gof{j}.rmse);
            H.Table.Data(4*(j-1)+1,3:4) = {r2, rmse}; 
            
            % plot measured and estimated vaues
            H.EstMarker(j,1) = line(H.AX,...
                H.Mod.DPFreq(j),H.Mod.DPMeasuredFrequency(j),...
                'markersize',10,'Marker','<','linestyle','none',...
                'Color',[0.1 0.6 0.3],'markerfacecolor',[0.1 0.6 0.3],...
                'linewidth', 1);
            H.EstMarker(j,2) = line(H.AX,...
                H.fa(j),H.Mod.DPEnergyAverage(j),...
                'markersize',10,'Marker','>','linestyle','none',...
                'Color',[0.6 0.3 0.1],'markerfacecolor',[0.6 0.3 0.1],...
                'linewidth', 1);
             H.EstMarker(j,3) = line(H.AX,...
                 H.fa(j),H.Mod.DPEstimGenCompAmp(j),...
                 'markersize',10,'Marker','v','linestyle','none',...
                 'Color',[0.3 0.1 0.6],'markerfacecolor',[0.3 0.1 0.6],...
                 'linewidth', 1);
             H.Legend = legend(H.AX,[H.L H.M(j,:) H.EstMarker(j,:)],...
                {'dpoae','noise','Model','Model Input',...
                'Measured Value','Energy Average','Model Estimate (A)'},...
                'Location','northwest');
            end
        end      
    end
    
%    set(H.Estimates.Children,'Enable','on')
    setappdata(H.mainfig,'H',H)
end

function EstimateSelection(hObject,EvntData)
H = getappdata(hObject.Parent,'H');
OldSelec = EvntData.OldValue;
NewSelec = EvntData.NewValue;


switch NewSelec
    case H.Est(1)
        H.EstMarker(1) = line(H.AX,H.Mod.DPFreq,H.Mod.DPMeasuredFrequency,'markersize',10,...
            'Marker','<','linestyle','none',...
            'Color',[0.1 0.6 0.3],'markerfacecolor',[0.1 0.6 0.3],'linewidth', 1);
        
    case H.Est(2)
        H.EstMarker(2) = line(H.AX,H.fa,H.Mod.DPEnergyAverage,'markersize',10,...
            'Marker','>','linestyle','none',...
            'Color',[0.6 0.3 0.1],'markerfacecolor',[0.6 0.3 0.1],'linewidth', 1);
            
    case H.Est(3)
        H.EstMarker(3) = line(H.AX,H.fa,H.Mod.DPEstimGenCompAmp,'markersize',10,...
            'Marker','v','linestyle','none',...
            'Color',[0.3 0.1 0.6],'markerfacecolor',[0.3 0.1 0.6],'linewidth', 1);
                
    case H.Est(4)
        delete(H.EstMarker);
        delete(H.M);
        set(H.Estimates.Children,'Enable','off')
end
setappdata(H.mainfig,'H',H)
end

%%
% FIG2 = figure(500);
% FIG2.PaperType = 'A4';
% FIG2.PaperOrientation = 'portrait';
% FIG2.PaperPosition = [0.25 0.25 7.768 11.193];
% AX1 = axes('parent',FIG2,'nextplot','add','ylim',[-20 20],...
%         'xlim',[850 6500],'xscale','log','box','on',...
%         'xtick',[1000 1500 2000 3000 4000 5000],'Xminortick','off');
% errorbar(AX1,fa,mean(DPEnergyAverage,2),std(DPEnergyAverage,[],2),'Color','b','linewidth', 2)
% errorbar(AX1,fa,mean(DPMeasuredFrequency,2),std(DPMeasuredFrequency,[],2),'Color','r','linewidth', 2)        
% errorbar(AX1,fa,mean(DPEstimGenCompAmp,2),std(DPEstimGenCompAmp,[],2),'Color',[0.3 1 0.3],'linewidth', 2)        
% K(1) = line(AX1,fa,mean(DPEnergyAverage,2),'marker','o','markersize',8,...
%             'color','b','markerfacecolor','b','linewidth', 1);
% K(2) = line(AX1,fa,mean(DPMeasuredFrequency,2),'marker','s','markersize',8,...
%             'color','r','markerfacecolor','r','linewidth', 1);
% K(3) = line(AX1,fa,mean(DPEstimGenCompAmp,2),'marker','d','markersize',8,...
%             'color',[0 0.8 0],'markerfacecolor',[0 0.8 0],'linewidth', 1);
% 
% legend(AX1,K,{'Energy Average', 'Measured', 'Estimated'},'location','best')
% 
% ylabel(AX1,'DPOAE amplitude [dB SPL]');
% xlabel(AX1,'Frequency [Hz]');   
% print(FIG2,'-dpdf', 'MeasVsModel.pdf');


