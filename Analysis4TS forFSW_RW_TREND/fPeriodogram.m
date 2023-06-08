function [ vPerio, vomega ] = fPeriodogram(vy)
% Periodogram via FFT
cn = length(vy); % number of observations
vomega = 2*pi*(0:cn-1)/cn; % fourier frequencies
vPerio = (abs(fft(vy-mean(vy))).^2)/(cn * 2 * pi) ; 
plot(vomega(2:end),vPerio(2:end),'b')
end


