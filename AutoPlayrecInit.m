function AutoPlayrecInit(fs)

%% Platform independent initialization of Playrec
% AUTOPLAYRECINIT(FS) automatically initialise playrec, if no
% input arguments is given, the default sampling frequency is 48kHz.

% change the following path to your Playrec path
%addpath /home/rop/installs/playrec_2_1_0/

if nargin<1; fs = 48e3; end

if playrec('isInitialised') == 0

    Devices = size(playrec('getDevices'));
    H = playrec('getDevices');

%    SC = 'Built-in Output';
     SC = 'ASIO Fireface USB';
%SC = 'QUAD-CAPTURE';
% %    SC = 'HDA Int';
H = playrec('getDevices');
if playrec('isInitialised') == 0

    for i = 1:1:size(H,2)
       if  strcmp(H(i).name,SC) 
           OutDevID = i-1;
           InDevID = i-1;
           break;
        end
    end
%    if  exist('OutDevID') == 1
        playrec('init',fs,OutDevID,InDevID,2,2,480,10e-3,10e-3)
%         fprintf('\nEdirol Sound Card (SC) is initialized\nIf "MON SW" lights, push it down\n\n')
%    % else
%         fprintf('\nEdirol Sound Card (SC) is not connected to the computer\nor push "MON SW" on SC and try again\nor restart Matlab and try again\n\n')
%     end

else
    fprintf('\nSound Card is already initialized\n\n')
end
end



