% Alexander Hay
% Homework 6

fprintf('\n');
fprintf('Alexander Hay\n');
fprintf('NUIN 408\n');
fprintf('Homework 6\n');

% α σ μ ≠ 
%% Problem 1
% ***********************************************
% 1a

fprintf('\n');
fprintf('Problem 1a *****************************\n');
fprintf('\n');

fprintf('Signal 1 matches spectrum b\n');
fprintf('Signal 2 matches spectrum c\n');
fprintf('Signal 3 matches spectrum a\n');
fprintf('\n');
fprintf('yoof\n');

%% Problem 2
% ***********************************************
% 2a

fprintf('\n');
fprintf('Problem 2a *****************************\n');
fprintf('\n');

load('HW6_prob2.mat');

threshold = 200;
fs = 10000;
timeaxis = 0:1/fs:(length(spikeTrace)-1)/fs;

[peak_sizes, peakPoints] = findpeaks(spikeTrace, 'MinPeakHeight', threshold);

figure_2a = figure;
hold('on');
plot(timeaxis, spikeTrace);
plot(timeaxis(peakPoints), peak_sizes, 'rx');
title('2a - Simple Treshold');
xlabel('Time (s)');
ylabel('Excitement (pA)');
hold off;

fprintf('see figure 2a\n');

%% **********************************************
% 2b

fprintf('\n');
fprintf('Problem 2b *****************************\n');
fprintf('\n');

fs = 10000/5;
y = downsample(spikeTrace,5);
timeaxis = 0:1/fs:(length(y)-1)/fs;
[peak_sizes, peakPoints] = findpeaks(y, 'MinPeakHeight', threshold);

figure_2b = figure;
hold on;
plot(timeaxis,y);
plot(timeaxis(peakPoints), peak_sizes, 'rx');
title('2b - Simple Treshold - Downsample');
xlabel('Time (s)');
ylabel('Excitement (pA)');
hold off;

fprintf('see figure 2b\n');
fprintf('\n');
fprintf('tweaks aside so the code runs (array size) it seems to have worked still.\n');

%% **********************************************
% 2c

fprintf('\n');
fprintf('Problem 2c *****************************\n');
fprintf('\n');

y = resample(spikeTrace,1,5);
timeaxis = 0:1/fs:(length(y)-1)/fs;
[peak_sizes, peakPoints] = findpeaks(y, 'MinPeakHeight', threshold);

figure_2c = figure;
hold on;
plot(timeaxis,y);
plot(timeaxis(peakPoints), peak_sizes, 'rx');
title('2c - Simple Treshold -Resample');
xlabel('Time (s)');
ylabel('Excitement (pA)');
hold off;

fprintf('see figure 2c\n');
fprintf('\n');
fprintf('Downsample takes every nth sample. Resample interpolates the data and resamples at a new rate.\n');

%% Problem 3
% ***********************************************
% 3a

fprintf('\n');
fprintf('Problem 3a *****************************\n');
fprintf('\n');

load('HW6_prob3.mat');
load('filters.mat');

% Matlab's FFT documentation helped with this problem

fs = 10000;                             % Sampling frequency
T = 1/fs;                               % Period
L = length(CC_data);                    % Signal length
t = 0:1/fs:(length(CC_data)-1)/fs;      % Time axis

Y = fft(CC_data_corrupted_A);           % FFT the signal
P2 = abs(Y/L);                          % double-sided spectrum
P1 = P2(1:L/2+1);                       % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure_3a = figure;
hold on;
subplot(1,2,2,'Parent',figure_3a);
plot(f,P1);
title('Figure 3a - FFT of corrupted data A');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
 
Y = fft(CC_data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(1,2,1,'Parent',figure_3a);
plot(f,P1);
title('Figure 3a - FFT of normal data');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');

fprintf('see figure 3a\n');
fprintf('\n');
fprintf('Using FFT its clear that theres an 80Hz signal\n');

%% ***********************************************
% 3b

fprintf('\n');
fprintf('Problem 3b *****************************\n');
fprintf('\n');

load('HW6_prob3.mat');

fs = 10000;                             % Sampling frequency
T = 1/fs;                               % Period
L = length(CC_data);                    % Signal length
t = 0:1/fs:(length(CC_data)-1)/fs;      % Time axis

Y = fft(CC_data_corrupted_B);           % FFT the signal
P2 = abs(Y/L);                          % double-sided spectrum
P1 = P2(1:L/2+1);                       % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

figure_3b = figure;
hold on;
subplot(1,2,2,'Parent',figure_3b);
plot(f,P1);
title('Figure 3a - FFT of corrupted data B');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
 
Y = fft(CC_data);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

subplot(1,2,1,'Parent',figure_3b);
plot(f,P1);
title('Figure 3b - FFT of normal data');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');

fprintf('see figure 3b');
fprintf('\n');
fprintf('Using FFT its clear that theres an 2kHz signal\n');

%% ***********************************************
% 3c

fprintf('\n');
fprintf('Problem 3c *****************************\n');
fprintf('\n');


% Filter parameters:
% Highpass

% FIR: Equiripple

% Minimum order

% Density factor: 20

% Units: Hz
% FS: 10000
% Fstop: 75
% Fpass: 85

% Units: dB
% Astop: 60
% Apass: 1

CC_data_filt_A = filtfilt(hi80.Numerator,1,CC_data_corrupted_A);

figure_3c_1 = figure;
plot(CC_data,'k','DisplayName', 'Original Data');
hold on;
plot(CC_data_corrupted_A,'r','DisplayName','Corrupted Data');
plot((-60+CC_data_filt_A),'b','DisplayName','Filtered Data');
title('Figure 3c - Filtered Corrupted Data A');
xlabel('Time (s)');
ylabel('Excitement (pA)');
legend;
hold off;

% Filter parameters:
% Lowpass

% FIR: Equiripple

% Minimum order

% Density factor: 20

% Units: Hz
% FS: 10000
% Fstop: 1950
% Fpass: 2050

% Units: dB
% Astop: 80
% Apass: 1

CC_data_filt_B = filtfilt(lo2k.Numerator,1,CC_data_corrupted_B);

figure_3c_2 = figure;
plot(CC_data_corrupted_B,'r','DisplayName','Corrupted Data');
hold on;
plot((-6.5+CC_data_filt_B),'b','DisplayName','Filtered Data');
plot(CC_data,'k','DisplayName', 'Original Data');
title('Figure 3c - Filtered Corrupted Data B');
xlabel('Time (s)');
ylabel('Excitement (pA)');
legend;
hold off;

fprintf('see figures 3c');
fprintf('\n');
fprintf('A lowpass filter for data A with limited success. A highpass filter was used for B.\n');
fprintf('A bandpass filter could be used for data B as well\n');

%% ***********************************************
% 3d

fprintf('\n');
fprintf('Problem 3d *****************************\n');
fprintf('\n');

[freqAxis, power] = powerSpectrum(CC_data, fs);
[freqAxis_corrupt_A, power_corrupt_A] = powerSpectrum(CC_data_corrupted_A, fs);
[freqAxis_corrupt_B, power_corrupt_B] = powerSpectrum(CC_data_corrupted_B, fs);
[freqAxis_A, power_A] = powerSpectrum(CC_data_filt_A, fs);
[freqAxis_B, power_B] = powerSpectrum(CC_data_filt_B, fs);

figure_3d = figure;
subplot(2,2,1,'Parent',figure_3d);
plot(CC_data_corrupted_A,'r','DisplayName','Corrupted Data');
hold on;
plot((-60+CC_data_filt_A),'b','DisplayName','Filtered Data');
plot(CC_data,'k','DisplayName','Original Data');
title('Data A - Time Domain');
xlabel('Time (s)');
ylabel('Excitement (pA)');
legend;

subplot(2,2,2,'Parent',figure_3d);
plot(CC_data_corrupted_B,'r','DisplayName','Corrupted Data');
hold on;
plot((-6.5+CC_data_filt_B),'b','DisplayName','Filtered Data');
plot(CC_data,'k','DisplayName','Original Data');
title('Data B - Time Domain');
xlabel('Time (s)');
ylabel('Excitement (pA)');
legend;

subplot(2,2,3,'Parent',figure_3d);
loglog(freqAxis_corrupt_A, power_corrupt_A,'r','DisplayName','Corrupted Data');
hold on;
loglog(freqAxis, power,'k','DisplayName','Original Data');
loglog(freqAxis_A, power_A,'b','DisplayName','Filtered Data');
title('Data A - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
legend;

subplot(2,2,4,'Parent',figure_3d);
loglog(freqAxis_corrupt_B, power_corrupt_B,'r','DisplayName','Corrupted Data');
hold on;
loglog(freqAxis, power,'k','DisplayName','Original Data');
loglog(freqAxis_B, power_B,'b','DisplayName','Filtered Data');
title('Data B - Frequency Domain');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
legend;

%% ***********************************************
% 3e

fprintf('\n');
fprintf('Problem 3d *****************************\n');
fprintf('\n');

% Attenuation = 10*log10(Pout./Pin)

att_a = 10*log10(power_corrupt_A./power_A);
att_b = 10*log10(power_corrupt_B./power_B);

figure_3e = figure;
subplot(1,2,1,'Parent',figure_3e);
plot(att_a);
hold on;
title('Figure 3e - Filter A Attenuation');
xlabel('Frequency (Hz)');
ylabel('Loss');

subplot(1,2,2,'Parent',figure_3e);
plot(att_b);
title('Figure 3e - Filter A Attenuation');
xlabel('Frequency (Hz)');
ylabel('Loss');
hold off;

fprintf('see figure 3e\n');
%% Problem 4
% ***********************************************
% 4a

fprintf('\n');
fprintf('Problem 4a *****************************\n');
fprintf('\n');

load('HW6_prob4.mat');


fs = 1000;
repeated_mean = mean(data_repeatedSeed);

[freqAxis_repeated, power_repeated] = powerSpectrum(repeated_mean, fs);
[freqAxis_nonRepeated, power_nonRepeated] = powerSpectrum(data_nonRepeatedSeed, fs);

% 

%% Problem 5
% ***********************************************
% 5a

fprintf('\n');
fprintf('Problem 5a *****************************\n');
fprintf('\n');

fs = 1000;            % Sampling frequency                    
T = 1/fs;             % Sampling period       
l = 2000;             % Length of signal
timeAxis = (0:l-1)*T;        % Time vector

signal_1 = sin(2*pi*10*(timeAxis));         % 10Hz sine wave, amplitude 1
signal_2 = .5*sin(2*pi*100*timeAxis);       % 100Hz sine wave, amplitude 0.5

filter_1 = movmean(signal_2,T);

figure_5a = figure;
subplot(3,1,1,'Parent',figure_5a);
plot(signal_1,'b','DisplayName','10Hz sine');
hold on;
plot(sys1(signal_1,fs),'r','DisplayName','sys1');
title('Function 1 vs 10Hz Sine Wave');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

subplot(3,1,2,'Parent',figure_5a);
plot(signal_2,'b','DisplayName','100Hz sine');
hold on;
plot(sys2(signal_2,fs),'r','DisplayName','sys2');
title('Function 2 vs 1kHz Sine Wave');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

subplot(3,1,3,'Parent',figure_5a);
plot(signal_1+0.25*randn(size(timeAxis)),'b','DisplayName','Noisy 10Hz sine');
hold on;
plot(sys3(signal_1,fs),'r','DisplayName','sys3');
title('Function 3 vs 10Hz Corrupted Sine Wave');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

fprintf('see figure 5a\n');
fprintf('\n');
fprintf('sys1 appears to be shifted by period T, but the signal starts curved rather than immediately upward.\n');
fprintf('sys2 is nonlinear; converges at time T.\n');
fprintf('sys3 appears to be the signal corrupted by some gaussian noise (demonstrated with noisy signal).\n');

%% **********************************************
% 5b

fprintf('\n');
fprintf('Problem 5b *****************************\n');
fprintf('\n');

s1 = sys1(signal_1,fs);
[freqAxis_1, power_1] = powerSpectrum(signal_1, fs);
[freqAxis_sys1, power_sys1] = powerSpectrum(s1, fs);

s2 = sys2(signal_2,fs);
[freqAxis_2, power_2] = powerSpectrum(signal_2, fs);
[freqAxis_sys2, power_sys2] = powerSpectrum(s2, fs);

s3 = sys3(signal_1,fs);
[freqAxis_3, power_3] = powerSpectrum(signal_1, fs);
[freqAxis_sys3, power_sys3] = powerSpectrum(s3, fs);

figure_5b = figure;
subplot(1,3,1,'Parent',figure_5b);
loglog(freqAxis_1, power_1,'b','DisplayName','10Hz sine');
hold on;
loglog(freqAxis_sys1, (10e-30)*power_sys1,'r','DisplayName','sys1');
title('5b - Function 1 Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');
legend;
hold off;

subplot(1,3,2,'Parent',figure_5b);
loglog(freqAxis_2, power_2,'b','DisplayName','100Hz sine');
hold on;
loglog(freqAxis_sys2, (10e-26)*power_sys2,'r','DisplayName','sys2');
title('5b - Function 2 Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');
legend;
hold off;

subplot(1,3,3,'Parent',figure_5b);
loglog(freqAxis_3, power_3,'b','DisplayName','100Hz sine');
hold on;
loglog(freqAxis_sys3, (10e-29)*power_sys3,'r','DisplayName','sys3');
title('5b - Function 3 Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');
legend;
hold off;

fprintf('Even though it phase shifts, sys1 creates higher resonant frequencies.\n');
fprintf('sys2 creates resonant frequencies around the input signal.\n');
fprintf('sys3 is a gaussian filter\n');

%% **********************************************
% 5c

fprintf('\n');
fprintf('Problem 5c *****************************\n');
fprintf('\n');

[freqAxis2, power2] = powerSpectrum(signal_2, fs);
% s2 = fft(sys2(signal_2,fs));

fs = 1000;                              % Sampling frequency                    
T = 1/fs;                               % Sampling period       
l = length(signal_2);                   % Length of signal
t = 0:1/fs:(length(signal_2)-1)/fs;     % Time axis

Y = fft(sys2(signal_2,fs));             % FFT the signal
P2 = abs(Y/l);                          % double-sided spectrum
P1 = P2(1:l/2+1);                       % single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(l/2))/l;

figure_5c1 = figure;
plot(signal_1,'b','DisplayName','Signal 1');
hold on;
plot(sys1(signal_1,fs),'r','DisplayName','sys1');
title('Function 1 vs 10Hz Sine Wave');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

figure_5c2 = figure;
subplot(1,3,1,'Parent',figure_5c2);
plot(signal_2,'b','DisplayName','Signal 1');
hold on;
plot(sys2(signal_2,fs),'r','DisplayName','sys1');
title('5c - Function 2 Time Domain');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

subplot(1,3,2,'Parent',figure_5c2);
loglog(freqAxis2, power2,'b','DisplayName','signal 2 Data');
hold on;
[freqAxis2, power2] = powerSpectrum(sys2(signal_2,fs), fs);
loglog(freqAxis2, (10e-25)*power2,'r','DisplayName','sys2 Data');
title('5c - Function 2 Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');
legend;
hold off;

subplot(1,3,3,'Parent',figure_5c2);
hold on;
plot(f,P1);
title('Figure 5c - FFT Function 2');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
legend;
hold off;

figure_5c3 = figure;
plot(signal_1+0.25*randn(size(timeAxis)),'b','DisplayName','Noisy Signal 1');
hold on;
plot(sys3(signal_1,fs),'r','DisplayName','sys3');
title('Function 3 vs 10Hz Corrupted Sine Wave');
xlabel('Time (ms)');
ylabel('Amplitude');
legend;
hold off;

fprintf('sys1 appears to be a phase shift, a linear function.\n');
fprintf('\n');
fprintf('The time it takes the signal from s2 to converge on its input signal is 1 period.\n');
fprintf('(visible in the time domain in figure 5c)\n');
fprintf('\n');
fprintf('sys3 is a gaussian filter, demonstrated by an input signal corrupted by gaussian noise.\n');