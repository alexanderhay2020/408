%% load data
load('noisySpikeData.mat')
Fs = 2E4; %sampling freq.
sampleInterval = 1./Fs;
Npoints = length(dataTrace);

timeAxis = [0:sampleInterval:sampleInterval*(Npoints-1)];

%% Plot trace and its power spectrum
[freqAxis, dataTrace_power] = powerSpectrum(dataTrace, Fs);

figure(1);
subplot(1,2,1);
plot(timeAxis, dataTrace, 'k');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Time domain');

subplot(1,2,2);
loglog(freqAxis, dataTrace_power, 'k');
xlabel('Frequency (Hz)');
ylabel('Power (uV^2/Hz)');
title('Frequency domain');

%% Try smoothing
span = 10;
data_smooth_span10 = smooth(dataTrace,span);
span = 30;
data_smooth_span30 = smooth(dataTrace,span);

figure(2);
plot(timeAxis, dataTrace, 'k');
hold('on');
plot(timeAxis, data_smooth_span10, 'r', 'linewidth', 2);
plot(timeAxis, data_smooth_span30, 'b', 'linewidth', 2);
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Time domain');
legend({'Original', 'Smoothed_ span 10', 'Smoothed_ span 30'});
hold('off');

%% Try high pass filtering
load('HP250.mat')
data_highPass_filtered = filtfilt(HighPass_250Hz.Numerator, 1, dataTrace);

figure(3);
plot(timeAxis, dataTrace, 'k');
hold('on');
plot(timeAxis, data_highPass_filtered, 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Time domain');
legend({'Original', 'High pass filtered'});
hold('off');

% %% Set a threshold and use it to find putative spikes
% threshold = 400;
% [peak_sizes, peakPoints] = findpeaks(-dataTrace, 'MinPeakHeight', threshold);
% 
% Npeaks = length(peakPoints);
% 
% window = -50:50;
% window_length = length(window);
% 
% spikeWaveformMatrix = zeros(Npeaks, window_length);
% for i=1:Npeaks
%     spikeWaveformMatrix(i,:) = dataTrace(peakPoints(i)+window);
% end
% 
% timeAxis_waveform = 1E3*window*sampleInterval;
% 
% figure(4);
% for i=1:Npeaks
%     plot(timeAxis_waveform, spikeWaveformMatrix(i,:), 'k');
%     pause;
% end

%% Back to the original data... let's set a threshold to find large events

figure(4);
plot(timeAxis, dataTrace, 'k');
xlabel('Time (s)');
ylabel('Voltage (uV)');


%% find peaks
%first flip the data over so we can look for positive peaks

threshold = 600;

figure(4);
plot(timeAxis, -dataTrace, 'k');
hold('on');
drawline('Position', [0 threshold; timeAxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Inverted data');
hold('off');

%% Find peaks
[peak_sizes, peakPoints] = findpeaks(-dataTrace, 'MinPeakHeight', threshold);

%% Plot them on the data
figure(5);
hold('on');
plot(timeAxis, -dataTrace, 'k');
plot(timeAxis(peakPoints), peak_sizes, 'rx');
drawline('Position', [0 threshold; timeAxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Inverted data');
hold('off');

%% Let's make a matrix and average the data around each of the peaks
% 50 points on either side of it
Npeaks = length(peakPoints);

window = -50:50;
window_length = length(window);

spikeWaveformMatrix = zeros(Npeaks, window_length);
for i=1:Npeaks
    spikeWaveformMatrix(i,:) = dataTrace(peakPoints(i)+window);
end

%% Get the mean and s.e.m. of the spikeWaveformMatrix and plot them
spikeWaveform_mean = mean(spikeWaveformMatrix, 1);
spikeWaveform_sem = std(spikeWaveformMatrix, [], 1) ./ sqrt(Npeaks);

timeAxis_waveform = 1E3*window*sampleInterval;

figure(6);
errorbar(timeAxis_waveform, spikeWaveform_mean, spikeWaveform_sem);
xlabel('Time from peak (ms)');
ylabel('Average voltage (uV)');
title('Average spike waveform for large peaks');

%% Now use this mean spike waveform as a filter

spikeWaveform_Filter = spikeWaveform_mean./sum(spikeWaveform_mean);

dataTrace_matchFiltered = filtfilt(spikeWaveform_Filter, 1, -dataTrace);
figure(7);
subplot(1,2,1);
plot(timeAxis, dataTrace, 'k');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('Original');
subplot(1,2,2);
plot(timeAxis, dataTrace_matchFiltered, 'k');
xlabel('Time (s)');
ylabel('Filtered data (arb. units');
title('Matched filtered');


%% Power spectrum of matched filtered data
[freqAxis, dataTrace_matchFiltered_power] = powerSpectrum(dataTrace_matchFiltered, Fs);

figure(8);
subplot(1,2,1);
plot(timeAxis, dataTrace_matchFiltered, 'k');
xlabel('Time (s)');
title('Time domain');

subplot(1,2,2);
loglog(freqAxis, dataTrace_matchFiltered_power, 'k');
xlabel('Frequency (Hz)');
title('Frequency domain');

%% Find peaks on the match filtered data
threshold = 7.5E3;
[peak_sizes, peakPoints] = findpeaks(dataTrace_matchFiltered, 'MinPeakHeight', threshold);

%% Collect them
Npeaks = length(peakPoints);

window = -50:50;
window_length = length(window);

spikeWaveformMatrix = zeros(Npeaks, window_length);
for i=1:Npeaks
    spikeWaveformMatrix(i,:) = dataTrace(peakPoints(i)+window);
end

%% Plot them
figure(9);
imagesc(spikeWaveformMatrix);
xlabel('Time from peak (ms)');
ylabel('Wavefornm number');
xticks(0:10:100)
xticklabels(timeAxis_waveform(1:10:end))
title('Identified spikes');
cmap = [0 0 1; 1 1 1; 1 0 0];
cmap_interp = interp1(-1:1:1, cmap, -1:.01:1);
colormap(cmap_interp);
caxis([-600 600]);
c = colorbar;
c.Label.String = 'Amplitude (uV)';

% figure(8);
% for i=1:Npeaks
%     plot(timeAxis_waveform, spikeWaveformMatrix(i,:), 'k');
%     pause;
% end






