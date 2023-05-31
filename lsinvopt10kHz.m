function [out] = lsinvopt10kHz(in)

% Inverse filter coefficients of Loudspeaker 1 with applied lowpass filter
% at 10 kHz
% ls1_inv_10kHz = [0.00171818136500944,0.0114197763665906,0.0416384442068358,...
%                  0.103737102884307,0.191317598569703,0.266310362869000,...
%                  0.272096528780350,0.178463626576396,0.0273161031811513,...
%                  -0.0770242282735591,-0.0477594924066765,0.0928645717941258,...
%                  0.210405188370008,0.170799769240083,-0.0241459675542668,...
%                  -0.221786591632874,-0.256837350859181,-0.113962091150319,...
%                  0.0571137528274939,0.0942392850291646,-0.0157069780435407,...
%                  -0.138287187503158,-0.143870454996370,-0.0376124957970973,...
%                  0.0698145423648839,0.0954514466140432,0.0599579273753730,...
%                  0.0329336817943382,0.0395615345805599,0.0469236642633420,...
%                  0.0289906738481464,0.00677977737412801,0.0119082604541982,...
%                  0.0344984841493992,0.0357233845367071,0.00537244697480024,...
%                  -0.0206045551187601,-0.0102791377684775,0.0184667486628305,...
%                  0.0224454391379753,-0.00577090588489761,-0.0273213730819711,...
%                  -0.0118959199037032,0.0194554280436481,0.0241222713076020,...
%                  -0.00243108084987508,-0.0211825837339353,-0.00507493241011300,...
%                  0.0236005809999600,0.0256695260211960,-1.66478292562625e-05];
ls1_inv_10kHz = [0.001718181365009440039401, 0.005500590845842727970305, ...
    0.011012914943324937644409, 0.012195171417047200035366, ...
    0.005580602604862153898191, -0.008013298798296798752805, ...
    -0.016991700847578977695562, -0.012426228893253177254508, ...
    0.005065186699650087165381, 0.018404353222247954446900, ...
    0.014969658901089977209864, -0.003931171193244258066368, ...
    -0.017329830138660049648758, -0.012499973300841401430139, ...
    0.005279124466795279135845, 0.014117696353806892528571,...
    0.005618718907974863674415, -0.008936692613696258294387, ...
    -0.010083010483518702396499, 0.002394552387258963924155,...
    0.012216329486471253634727, 0.005802828324430232137532, ...
    -0.007909785988360697878141, -0.011439758228434319595190,...
    -0.000035401337465346474631, 0.010658103408513397089563, ...
    0.007186986441465796453254, -0.005490534973999555527768, ...
    -0.010066163819160183628965, -0.001337182440680927640597, ...
    0.008534531319412931582580, 0.006102833003204326776936, ...
    -0.004462612394483264979805, -0.008866431577861905166671, ...
    -0.001400943605660446254768, 0.007117343105861254591582,...
    0.005599767669254647646515, -0.003270585420619757238692,...
    -0.007084955934624826735801, -0.001542549473679133940907,...
    0.005143063120944644725507, 0.003987727213311449171729,...
    -0.002265821409167293690873, -0.004602398085276322194093,...
    -0.000390362059353477320288, 0.003550119821060727237638,...
    0.001751824881594757398265, -0.002507170737615221477873,...
    -0.002647754177055239386412, 0.001271412338946530974637,...
    0.003165782491432467867648];




out = filter(ls1_inv_10kHz,1,in);

end
             
%% Data to this function is generated by the code below            
% clc; clear all; close all;
% 
% addpath ../../../Misc/Anders/lab/log;
% 
% fs = 48000;
% 
% load('mls-1rep-4avg-16-480-coupler-LS1L-LS1LMIC-h_aver1isCMic-h_aver2isPMic.mat','h_aver1');
% ls1 = h_aver1;
% ls1 = 23.7411.*ls1;
% 
% nfft = length(ls1);
% 
% order = 50;
% [a1, gain] = lpc(double(ls1),order);
% a1 = 1.5299.*a1;
% 
% [B A] = ellip(12,0.1,100,[10e3]/(fs/2),'low');
% a1 = filter(B,A,a1);

% 
% fs = 48000;
% 
% Measurements made with Impulse response measurements app. 
% ls1 = impresp.Amplitude;
% ls1 = 23.7411.*ls1;
% 
% nfft = length(ls1);
% 
% order = 50;
% [a1, gain] = lpc(double(ls1),order);
% a1 = 1.5299.*a1;
% 
% [B A] = ellip(12,0.1,100,[10e3]/(fs/2),'low');
% a1 = filter(B,A,a1);
% 






