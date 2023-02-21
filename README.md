# OAEsystem
 Matlab code for OAE measurements
This code was developed in MATLAB and uses playrec http://www.playrec.co.uk/ for audio playback and recording. The file playrec.mexw64 is a mex file for matlab that implements the playrec functionality for windows 64bit systems.

The program is made for a Etymotic Reseach ER10C OAE probe using a USB sound card using ASIO drivers in windows. Minimal requirements for the sound card are sampling frequency of 48 kHz and full duplex, two channel out and one channel in. The calibration is made using sound pressure level an artificial ear IEC 711 coupler with the correct adaptor.
