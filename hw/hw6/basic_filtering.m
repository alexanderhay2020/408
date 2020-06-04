%% filters and filter designs

load('filter_examples.mat')

%% we are looking for spikes, they look like this
figure()
plot(d)

% these are the 2 datasets we have, each has some kind of noise in it
fs = 1000;
timeaxis = 0:1/fs:(length(A)-1)/fs;

figure()
plot(timeaxis,A)
xlabel('Seconds');
ylabel('Voltage (uV)');
title('Trace A');
figure()
plot(timeaxis,B)
xlabel('Seconds');
ylabel('Voltage (uV)');
title('Trace B');
%% find peaks by eye by drawing a horizontal line and looking when the peaks are above it

threshold = 2;

figure(4);
plot(timeaxis, A, 'k');
hold('on');
drawline('Position', [0 threshold; timeaxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
hold('off');

%% Matlab has a great function that works the same way
[peak_sizes, peakPoints] = findpeaks(A, 'MinPeakHeight', threshold);

%% lets plot our found peaks with red x's
figure();
hold('on');
plot(timeaxis, A, 'k');
plot(timeaxis(peakPoints), peak_sizes, 'rx');
% drawline('Position', [0 threshold; timeaxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('trace A');
hold('off');

%% remember we are looking for spikes that look like this
figure()
plot(d)

%% We don't know which of those peaks are spikes and which are noise, so how do we fix that?

cleanA = lowpass(A,200,1000)

figure()
hold on
plot(A(1:100), 'k')
plot(cleanA(1:100), 'r')
hold off

%% now lets find peaks and plot our new trace
threshold = 1;
[peak_sizes, peakPoints] = findpeaks(cleanA, 'MinPeakHeight', threshold);
figure();
hold('on');
plot(timeaxis, cleanA, 'k');
plot(timeaxis(peakPoints), peak_sizes, 'rx');
drawline('Position', [0 threshold; timeaxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('trace A');
hold('off');

%% What about for B?

plot(timeaxis,B)
xlabel('Seconds');
ylabel('Voltage (uV)');
title('Trace B');


%% B has a drift in it, we cant just draw a threshold line, so what do we do?

cleanB = highpass(B, 1, fs);

figure()
hold on
plot(B(1:100), 'k')
plot(cleanB(1:100), 'r')
hold off

%% plot clean B
plot(timeaxis,cleanB)
xlabel('Seconds');
ylabel('Voltage (uV)');
title('Trace B');

%% now we can find our peaks and mark them as before
threshold = 1;
[peak_sizes, peakPoints] = findpeaks(cleanB, 'MinPeakHeight', threshold);
figure();
hold('on');
plot(timeaxis, cleanB, 'k');
plot(timeaxis(peakPoints), peak_sizes, 'rx');
drawline('Position', [0 threshold; timeaxis(end) threshold], 'Color', 'r');
xlabel('Time (s)');
ylabel('Voltage (uV)');
title('trace A');
hold('off');


%% filter designer
%when making an FIR filter, export it as either the coefficient (Num) or as
%the object (Hd)

%when making an IIR filter, export it as either the coefficients (SOS and
%G) or as the object (iirHd)

%% lets use the filter we created

lowpassedA = filtfilt(Num,1, A); %if FIR filter and exported coefficients
lowpassedA = filtfilt(Hd.Numerator,1,A); %if FIR filter and exported object

lowpassedA = filtfilt(SOS,G,A); %if IIR filter and exported coefficients
lowpassedA = filtfilt(iirHd.sosMatrix, iirHd.ScaleValues, A); % if IIR filter and exported object


figure()
hold on
plot(timeaxis, A, 'k')
plot(timeaxis, lowpassedA, 'r')
hold off

