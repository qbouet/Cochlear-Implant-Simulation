% TASK 1
% Load audio
[y,fs_original]=audioread('mbet.wav');
sound(y,fs_original);

% Calculate the duration of the audio file
duration = length(y) / fs_original; % Calculate the duration in seconds
t_original = (0:1/fs_original:duration)*1000; % Original time vector (484.4ms)

% Choose new sampling frequency and update the time vector
fs_new = 16000;
t_new = (0:1/fs_new:duration)*1000; % Desired time vector to match new sampling frequency (484.4ms)

% Resample at new sampling frequency (16kHz)
y1 = resample(y,fs_new,fs_original);

% Define center frequencies and bandwidths
center_frequencies = [250 375 500 625 750 937 1187 1437 1687 2000 2375 2812 3312 3875 4563 5375];
bandwidths = [124 126 124 126 124 250 250 250 250 376 374 500 500 626 750 874];

% Design low pass filter
Fc = 200; % Cutoff frequency (Hz)
N = 101; % Filter order (odd for symmetric FIR)

% Design the filter using the FIR filter design function
LPF = fir1(N, Fc/(fs_new/2), 'low');

% Apply a bandpass filter to each channel respectively
channels = [];
for i = 1:16
    % Get lower and upper frequency
    fL = center_frequencies(i) - 0.5*bandwidths(i);
    fH = center_frequencies(i) + 0.5*bandwidths(i);

    % Rectify signal and apply bandpass filter
    channel = abs(bandpass(y1,[fL,fH],fs_new));

    % Apply low pass filter
    channel = filter(LPF, 1, channel).';
    channels=[channels; channel];
end

% Normalise the signal
channels = channels./max(channels,[],"all");

% Map the amplitudes using the mapping function
C = 1;
T = 0;
A = (C - T)*channels + T;

% Sample every 4ms (equivalent to every 64th value for Fs of 16kHz)
channels2 = A(:,1:64:end);

% Select the 6 largest filter outputs every 4ms
channels3 = NaN(size(channels2));
for i = 1:length(channels2)
    % Sort the amplitudes
    [sorted_case, sorted_channel_indices] = sort(channels2(:,i), 'descend');
    % Select the highest 6 values from sorted amplitudes
    for k = 1:6
        channels3(sorted_channel_indices(k),i) = sorted_case(k);
    end
end

% Update the time vector
t_new2 = (0:64/(fs_new):duration)*1000; % (484.4ms)

% Plot electrodogram
figure;
tlo = tiledlayout("vertical", 'Padding', 'none', 'TileSpacing', 'none');
for i = 16:-1:2
    nexttile(tlo);
    plot(t_new2,channels3(i,:),'.','MarkerSize',15)
    set(gca,'XTick',[], 'YTick', [])
    ylabel(i)
    ylim([0,1]);
end
nexttile(tlo);
plot(t_new2,channels3(1,:),'.','MarkerSize',15)
set(gca, 'YTick', [])
ylabel(1)
ylim([0,1]);


% TASK 2
channels4 = [];
figure;

% Iterate over each spike
for i = 1:16
    % Use the filter's center frequency for the channel's waveform
    frequency = center_frequencies(i);

    % Interpolate amplitudes of channel to remove NaN values
    amplitude = fillmissing(channels3(i,:),'nearest'); % Choose 'nearest' method

    % Interpolate amplitudes again with time vector to achieve correct timing
    original_time = linspace(0, 484.4, length(amplitude)); % Initial time vector
    desired_time = linspace(0, 484.4, length(A)); % Desired time vector with correct length
    amplitude = interp1(original_time, amplitude, desired_time, 'pchip'); % Choose 'pchip' method

    % Generate waveform for channel
    sinusoid = amplitude.*sin(2*pi*[1:length(amplitude)]*frequency/fs_new);
    channels4(i,:) = sinusoid;
end

% Sum up generated waveforms
y2 = sum(channels4,1);
% Normalise the reconstructed signal
y2 = y2./max(y2);

% Plot and listen to reconstructed signal
plot(desired_time,y2)
hold on;
plot(t_new,y1)
sound(y2,fs_new);
title("Task 2 - Reconstructed Signal")
legend("reconstructed signal","original signal")


% TASK 3
channels5 = [];
figure;
% Define standard deviation for Gaussian function
std = 100/3;

% Iterate over each spike like in Task 2
for i = 1:16
    
    % Define shifts from -100 to 100 around center frequency
    shifts = -100:1:100;

    % Use the filter's center frequency for the channel's waveform
    center_frequency = center_frequencies(i);

    % Interpolate amplitudes of channel to remove NaN values
    amplitude = fillmissing(channels3(i,:),'nearest'); % Choose 'nearest' method

    % Interpolate amplitudes again with time vector to achieve correct timing
    original_time = linspace(0, 484.4, length(amplitude)); % Initial time vector
    desired_time = linspace(0, 484.4, length(A)); % Desired time vector with correct length
    amplitude = interp1(original_time, amplitude, desired_time, 'pchip'); % Choose 'pchip' method

    sinusoids = zeros(size(amplitude));
    for shift = shifts
        % Apply shift
        frequency = center_frequency + shift;
        mean = frequency;

        % Get factor and normalise using it by making the center frequency's amplitude "1"
        factor = (1/(std*sqrt(2*pi)))*exp((-0.5*(mean-shift-mean)^2)/std^2)/(1/(std*sqrt(2*pi)))*exp((-0.5*(mean-mean)^2)/std^2);
        % Adjust amplitude to the Gaussian function
        amplitude = factor*amplitude;

        % Generate waveform for channel
        sinusoid = amplitude.*sin(2*pi*[1:length(amplitude)]*frequency/fs_new);
        % Sum up sinusoids for the channel
        sinusoids = sinusoids + sinusoid;
    end
    channels5(i,:) = sinusoids;
end

% Sum up generated waveforms
y3 = sum(channels5,1);
% Normalise the reconstructed signal
y3 = y3./max(y3);

% Plot and listen to reconstructed signal
plot(desired_time,y3)
hold on;
plot(t_new,y1)
sound(y3,fs_new);
title("Task 3 - Reconstructed Signal")
legend("reconstructed signal","original signal")

