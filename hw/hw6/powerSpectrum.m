function [freq, power] = powerSpectrum(signal, sampleRate, bins)
L = length(signal);

if nargin>2
    freq = (0:sampleRate/bins:(sampleRate/2 - sampleRate/bins)); %set up frequency axis
    signal_omega = fft(signal, bins);
    power = signal_omega.*conj(signal_omega); %square output (removes imaginary part)
    power = power(1:round(bins/2))/bins; %take first side of power spectrum (symmetric for real signals)    
else
    freq = sampleRate*(0:(L/2))/L; %set up frequency axis
    signal_omega = fft(signal);
    power = signal_omega.*conj(signal_omega); %square output (removes imaginary part)
    power = power(1:round(L/2)+1)/L; %take first side of power spectrum (symmetric for real signals)
end
