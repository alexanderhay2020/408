%% load same data from last class and save key variables for each signal
load('powerSpectrumData');

% sample frequency and interval
sampleFreq = 1E4; %Hz
sampleInterval = 1/sampleFreq; %s

%got number of points
Npoints_sp = length(spikeData);
Npoints_CC = length(CC_data);
Npoints_VC = length(VC_data);

%make time axis (s)
timeAxis_sp = 0:sampleInterval:(Npoints_sp-1)*sampleInterval;
timeAxis_CC = 0:sampleInterval:(Npoints_CC-1)*sampleInterval;
timeAxis_VC = 0:sampleInterval:(Npoints_VC-1)*sampleInterval;

%% collect the power spectrum for each signal and plot it next to the signal
[freqAxis_sp, power_sp] = powerSpectrum(spikeData, sampleFreq);
[freqAxis_CC, power_CC] = powerSpectrum(CC_data, sampleFreq);
[freqAxis_VC, power_VC] = powerSpectrum(VC_data, sampleFreq);

figure(1);
subplot(3,2,1);
plot(timeAxis_sp, spikeData);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Cell-attached spike signal')
subplot(3,2,2);
loglog(freqAxis_sp, power_sp)
xlabel('Frequency (Hz)');
ylabel('Power (pA^2/Hz)');


subplot(3,2,3);
plot(timeAxis_CC, CC_data);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Current clamp signal')
subplot(3,2,4);
loglog(freqAxis_CC, power_CC)
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');

subplot(3,2,5);
plot(timeAxis_VC, VC_data);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Voltage clamp signal')
subplot(3,2,6);
loglog(freqAxis_VC, power_VC)
xlabel('Frequency (Hz)');
ylabel('Power (pA^2/Hz)');

%% Use a sliding window to smooth the VC data
span = 50; %number of points in moving average
VC_data_smoothed = smooth(VC_data, span);

figure(2);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_smoothed, 'r', 'linewidth', 2);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Voltage clamp signal')
legend({'Original', 'Smoothed'});

%% A moving average is the same as the convolution with a flat filter
F = ones(span,1) ./ span;
VC_data_conv_filter = conv(VC_data, F);

figure(3);
plot(timeAxis_VC, VC_data_conv_filter) 
%this returns an error because convolved signal is too long
length(timeAxis_VC)
length(VC_data_conv_filter)

%% plot the convolution
figure(3);
plot(VC_data_conv_filter)
title('Convolution of data and filter');

%% Remove the ends by calling conv with the 'same' argument
VC_data_conv_filter_same = conv(VC_data, F, 'same');
length(timeAxis_VC)
length(VC_data_conv_filter_same)

%% Now we can plot it on top of the data
figure(4);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_smoothed, 'r', 'linewidth', 2);
plot(timeAxis_VC, VC_data_conv_filter_same, 'g', 'linewidth', 1);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Voltage clamp signal')
legend({'Original', 'Matlab smooth function', 'Conv with flat filter'});

%% A longer span will filter more high frequencies out
VC_data_smoothed_span50 = VC_data_smoothed;

span = 200; %number of points in moving average
VC_data_smoothed_span200 = smooth(VC_data, span);

figure(5);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_smoothed_span50, 'r', 'linewidth', 1);
plot(timeAxis_VC, VC_data_smoothed_span200, 'g', 'linewidth', 1);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Voltage clamp signal')
legend({'Original', 'Smooth span 50', 'Smooth span 200'});

%% Let's look at these same traces in the frequency domain
[~, power_VC_span50] = powerSpectrum(VC_data_smoothed_span50, sampleFreq);
% ~ means don't bother returning this value
[~, power_VC_span200] = powerSpectrum(VC_data_smoothed_span200, sampleFreq);

figure(6);
subplot(1,2,1);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_smoothed_span50, 'r', 'linewidth', 1);
plot(timeAxis_VC, VC_data_smoothed_span200, 'g', 'linewidth', 1);
xlabel('Time (s)');
ylabel('Current (pA)');
title('Voltage clamp signal')
legend({'Original', 'Smooth span 50', 'Smooth span 200'});
hold('off');

subplot(1,2,2);
loglog(freqAxis_VC, power_VC)
hold('on');
loglog(freqAxis_VC, power_VC_span50, 'r')
loglog(freqAxis_VC, power_VC_span200, 'g')
xlabel('Frequency (Hz)');
ylabel('Power (pA^2/Hz)');
legend({'Original', 'Smooth span 50', 'Smooth span 200'});
hold('off');

%% Gaussian Filter
%The filter we convolve with the data does not have to be flat. 
%Let's try a gaussian
span = 50;
F = gausswin(span);
F = F./sum(F); %in order not to change the mean, the filter should sum to 1
figure(3);
plot(F);
xlabel('Data point');
ylabel('Filter coefficient');

%% Do the convolution and get its power spectrum
VC_data_conv_filter_gauss = conv(VC_data, F, 'same');
[~, power_VC_meanFilt] = powerSpectrum(VC_data_conv_filter_same, sampleFreq);
[~, power_VC_gaussFilt] = powerSpectrum(VC_data_conv_filter_gauss, sampleFreq);

%% plot the convolution
figure(4);
subplot(2,2,1);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_conv_filter_same, 'r', 'linewidth', 2);
plot(timeAxis_VC, VC_data_conv_filter_gauss, 'g', 'linewidth', 2);
legend({'Original', 'Mean filter', 'Gaussian filter'});
title('Full trace');
hold('off');

subplot(2,2,2);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_conv_filter_same, 'r', 'linewidth', 2);
plot(timeAxis_VC, VC_data_conv_filter_gauss, 'g', 'linewidth', 2);
legend({'Original', 'Mean filter', 'Gaussian filter'});
xlim([0.5, 0.7]);
title('Zoomed in');
hold('off');

subplot(2,2,3);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_conv_filter_same, 'r', 'linewidth', 2);
plot(timeAxis_VC, VC_data_conv_filter_gauss, 'g', 'linewidth', 2);
legend({'Original', 'Mean filter', 'Gaussian filter'});
xlim([0.5, 0.525]);
title('Zoomed in more');
hold('off');

subplot(2,2,4);
loglog(freqAxis_VC, power_VC);
hold('on');
loglog(freqAxis_VC, power_VC_meanFilt, 'r');
loglog(freqAxis_VC, power_VC_gaussFilt, 'g');
legend({'Original', 'Mean filter', 'Gaussian filter'});
xlabel('Freqency (Hz');
ylabel('Power (pA^2/Hz)');
title('Power spectrum');
hold('off');

%% CC data
%Now let's go to the current clamp data and try to filter it two different
%ways: one to highlight the spikes and one to highlight the subthreshold
%data and remove the spikes.
figure(7);
plot(timeAxis_CC, CC_data);
xlabel('Time (s)'); 
ylabel('Voltage (mV)');

%% Mean filter
span = 50;
CC_data_smoothed_mean = smooth(CC_data, span);
figure(8);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_smoothed_mean, 'r', 'linewidth', 2);
xlabel('Time (s)');
ylabel('Voltage (mV)');
hold('off');
legend({'Original', 'Mean filter'});

%% Median filter
span = 50;
CC_data_smoothed_median = smoothdata(CC_data, 'movmedian', span);

figure(9);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_smoothed_mean, 'r', 'linewidth', 2);
plot(timeAxis_CC, CC_data_smoothed_median, 'g', 'linewidth', 1);
xlabel('Time (s)');
ylabel('Voltage (mV)');
hold('off');
legend({'Original', 'Mean filter', 'Median filter'});

%% Look at the power spectrum of each signal
[~, power_CC_meanFilt] = powerSpectrum(CC_data_smoothed_mean, sampleFreq);
[~, power_CC_medianFilt] = powerSpectrum(CC_data_smoothed_median, sampleFreq);

figure(10);
subplot(1,2,1);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_smoothed_mean, 'r', 'linewidth', 2);
plot(timeAxis_CC, CC_data_smoothed_median, 'g', 'linewidth', 1);
xlabel('Time (s)');
ylabel('Voltage (mV)');
hold('off');
legend({'Original', 'Mean filter', 'Median filter'});

subplot(1,2,2);
loglog(freqAxis_CC, power_CC)
hold('on');
loglog(freqAxis_CC, power_CC_meanFilt, 'r')
loglog(freqAxis_CC, power_CC_medianFilt, 'g')
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');
legend({'Original', 'Mean filter', 'Median filter'});
hold('off');

%% Design a high pass filter to keep the spikes only
% use filterDesignApp
%F_highpass = highPass_100;
CC_data_highPass = filter(highPass_100Hz, CC_data);

%% Plot the result
figure(11);
plot(timeAxis_CC, CC_data_highPass);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('High pass filtered signal');

%% Use filtfilt to filter forwards and backwards to get rid of the edge artifact
CC_data_highPass = filtfilt(highPass_100Hz.Numerator, 1, CC_data);

%% plot it
figure(12);
plot(timeAxis_CC, CC_data_highPass);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('High pass filtered signal');

%% Design a Low Pass Filter to remove the spikes
%F_lowpass = lowPass_50;
CC_data_lowPass = filter(lowPass_50, CC_data);

%% plot it
figure(13);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_lowPass, 'r');
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Low pass filtered signal');
legend({'Original', 'Low pass filter'});

%% Use filtfilt
CC_data_lowPass = filtfilt(lowPass_50.Numerator, 1, CC_data);

%% plot it
figure(14);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_lowPass, 'r');
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Low pass filtered signal');
legend({'Original', 'Low pass filter'});

%% fix the offset
CC_data_lowPass = CC_data_lowPass + mean(CC_data) - mean(CC_data_lowPass);
figure(15);
plot(timeAxis_CC, CC_data);
hold('on');
plot(timeAxis_CC, CC_data_lowPass, 'r');
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Low pass filtered signal');
legend({'Original', 'Low pass filter'});

%% Plot each of the signals with its power spectrum
[~, power_CC_highpass] = powerSpectrum(CC_data_highPass, sampleFreq);
[~, power_CC_lowpass] = powerSpectrum(CC_data_lowPass, sampleFreq);

figure(16);
subplot(3,2,1);
plot(timeAxis_CC, CC_data);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Original');

subplot(3,2,2);
loglog(freqAxis_CC, power_CC);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');

subplot(3,2,3);
plot(timeAxis_CC, CC_data_highPass);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('High Pass Filtered');

subplot(3,2,4);
loglog(freqAxis_CC, power_CC_highpass);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');

subplot(3,2,5);
plot(timeAxis_CC, CC_data_lowPass);
xlabel('Time (s)');
ylabel('Voltage (mV)');
title('Low Pass Filtered');

subplot(3,2,6);
loglog(freqAxis_CC, power_CC_lowpass);
xlabel('Frequency (Hz)');
ylabel('Power (mV^2/Hz)');

%% Return to VC data and use a notch (band-stop) filter to remove the 60Hz noise
%F_notch = notch_60;
VC_data_notchFilt = filtfilt(notch_60.Numerator, 1, VC_data);

[~, power_VC_notchFilt] = powerSpectrum(VC_data_notchFilt, sampleFreq);
index = find(freqAxis_VC == 60)
power_VC(index)
power_VC_notchFilt(index)

%% Plot 
figure(17);
subplot(1,2,1);
plot(timeAxis_VC, VC_data);
hold('on');
plot(timeAxis_VC, VC_data_notchFilt, 'r');
xlabel('Time (s)');
ylabel('Voltage (mV)');
legend({'Original', 'Notch filtered'});
hold('off');

subplot(1,2,2);
loglog(freqAxis_VC, power_VC)
hold('on');
loglog(freqAxis_VC, power_VC_notchFilt, 'r')
xlabel('Frequency (Hz)');
ylabel('Power (pA^2/Hz)');
legend({'Original', 'Notch filtered'});
hold('off');
