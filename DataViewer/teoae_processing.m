function [repro RMS noiseRMS spectrum noise_spectrum mean_response noise A_av B_av] = teoae_processing(A,B)
    
A_av=mean(A,2);
B_av=mean(B,2);

repro=abs(corr(A_av,B_av))*100; % percentage

A_av=A_av'; B_av=B_av';
mean_response=(A_av+B_av)/2;
noise=(A_av-B_av)/2;

r=find(mean_response>0);
RMS=20*log10(sqrt(mean(mean_response(r(1):end).^2))/20e-6);
noiseRMS=20*log10(sqrt(mean(noise(r(1):end).^2))/20e-6);

nfft1=length(mean_response); 
spectrum=2*abs(fft(mean_response))/nfft1;
spectrum=20*log10(spectrum(1:nfft1/2)/20e-6);
noise_spectrum=2*abs(fft(noise,nfft1))/nfft1;
noise_spectrum=20*log10(noise_spectrum(1:nfft1/2)/20e-6);


end