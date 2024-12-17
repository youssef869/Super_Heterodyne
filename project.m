
clc; clear; close all;

% File paths for the input signals
file1 = 'Short_QuranPalestine.wav';
file2 = 'Short_FM9090.wav';
file3 = 'Short_BBCArabic2.wav';
% file4 = 'Short_RussianVoice.wav';
% file5 = 'Short_SkyNewsArabia.wav';

% Read the stereo signals
[stereo1, fs1] = audioread(file1);
[stereo2, fs2] = audioread(file2);
[stereo3, fs3] = audioread(file3);
% [stereo4, fs4] = audioread(file4);
% [stereo5, fs5] = audioread(file5);

% Convert stereo to mono by averaging the two channels
mono1 = stereo1(:,1) + stereo1(:,2);
mono2 = stereo2(:,1) + stereo2(:,2);
mono3 = stereo3(:,1) + stereo3(:,2);
% mono4 = stereo4(:,1) + stereo4(:,2);
% mono5 = stereo5(:,1) + stereo5(:,2);

% Zero-pad the shorter signal to match the length of the longer one
maxLength = max([length(mono1), length(mono2),length(mono3)]);
mono1 = [mono1; zeros(maxLength - length(mono1), 1)];
mono2 = [mono2; zeros(maxLength - length(mono2), 1)];
mono3 = [mono3; zeros(maxLength - length(mono3), 1)];
% mono4 = [mono4; zeros(maxLength - length(mono4), 1)];
% mono5 = [mono5; zeros(maxLength - length(mono5), 1)];

% Create time vectors for mono signals
t1_orig = (0:length(mono1)-1) / fs1; % Time vector for mono1
t2_orig = (0:length(mono2)-1) / fs2; % Time vector for mono2
t3_orig = (0:length(mono3)-1) / fs3; % Time vector for mono3
% t4_orig = (0:length(mono4)-1) / fs4; % Time vector for mono4
% t5_orig = (0:length(mono5)-1) / fs5; % Time vector for mono5

% % --- Plot Signal 1 BEFORE Upsampling: Time-Domain and Frequency-Domain ---
% figure;
% subplot(2,1,1);
% plot(t1_orig, mono1);
% title('Time-Domain Plot of Signal 1 (Before Upsampling)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % FFT for mono1
% N1 = length(mono1); % Length of signal
% Y1 = fft(mono1); % Compute FFT
% f1 = (-N1/2:N1/2-1) * (fs1 / N1); % Frequency axis for double-sided spectrum
% Y1_shifted = fftshift(Y1); % Shift the FFT for proper double-sided display
% 
% subplot(2,1,2);
% plot(f1, abs(Y1_shifted)); % Double-sided spectrum
% title('Spectrum of Signal 1 (Before Upsampling)');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;
% 
% % --- Plot Signal 2 BEFORE Upsampling: Time-Domain and Frequency-Domain ---
% figure;
% subplot(2,1,1);
% plot(t2_orig, mono2);
% title('Time-Domain Plot of Signal 2 (Before Upsampling)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % FFT for mono2
% N2 = length(mono2); % Length of signal
% Y2 = fft(mono2); % Compute FFT
% f2 = (-N2/2:N2/2-1) * (fs2 / N2); % Frequency axis for double-sided spectrum
% Y2_shifted = fftshift(Y2); % Shift the FFT for proper double-sided display
% 
% subplot(2,1,2);
% plot(f2, abs(Y2_shifted)); % Double-sided spectrum
% title('Spectrum of Signal 2 (Before Upsampling)');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;


% Define the upsampling factor
r = 10; % Increase sampling rate by 10 times

% Calculate new sampling frequency
new_fs = r * fs1; % New sampling frequency (44100 * 10 = 441000 Hz)

% Resample signals
mono1_resampled = interp(mono1, r); % Upsample mono1 by factor r
mono2_resampled = interp(mono2, r); % Upsample mono2 by factor r
mono3_resampled = interp(mono3, r); % Upsample mono3 by factor r
% mono4_resampled = interp(mono4, r); % Upsample mono4 by factor r
% mono5_resampled = interp(mono5, r); % Upsample mono5 by factor r

N_resampled = max([length(mono1_resampled), length(mono2_resampled), length(mono3_resampled)]);
% Create new time vectors for resampled signals
t_new = (0:N_resampled-1) / new_fs; % New time axis for mono
f_resampled = (-N_resampled/2:N_resampled/2-1) * (new_fs / N_resampled);


% --- Plot Signal 1 AFTER resampling: Time-Domain and Frequency-Domain ---
figure;
subplot(2,1,1);
plot(t_new, mono1_resampled);
title('Time-Domain Plot of Resampled Signal 1 (After resampling)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT for resampled mono1
N_resampled = length(mono1_resampled);
Y1_resampled = fft(mono1_resampled);
Y1_resampled_shifted = fftshift(Y1_resampled);

subplot(2,1,2);
plot(f_resampled, abs(Y1_resampled_shifted));
title('Spectrum of Resampled Signal 1 (After resampling)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% --- Plot Signal 2 AFTER resampling: Time-Domain and Frequency-Domain ---
figure;
subplot(2,1,1);
plot(t_new, mono2_resampled);
title('Time-Domain Plot of Resampled Signal 2 (After resampling)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT for resampled mono2
modulated2_freq = fft(mono2_resampled);
modulated2_freq_shifted = fftshift(modulated2_freq);

subplot(2,1,2);
plot(f_resampled, abs(modulated2_freq_shifted));
title('Spectrum of Resampled Signal 2 (After resampling)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% --- Plot Signal 3 AFTER resampling: Time-Domain and Frequency-Domain ---
figure;
subplot(2,1,1);
plot(t_new, mono3_resampled);
title('Time-Domain Plot of Resampled Signal 3 (After resampling)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;


modulated3_freq = fft(mono3_resampled);
modulated3_freq_shifted = fftshift(modulated3_freq);

subplot(2,1,2);
plot(f_resampled, abs(modulated3_freq_shifted));
title('Spectrum of Resampled Signal 3 (After resampling)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% % --- Plot Signal 4 AFTER resampling: Time-Domain and Frequency-Domain ---
% figure;
% subplot(2,1,1);
% plot(t_new, mono4_resampled);
% title('Time-Domain Plot of Resampled Signal 4 (After resampling)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% modulated4_freq = fft(mono4_resampled);
% modulated4_freq_shifted = fftshift(modulated4_freq);
% 
% subplot(2,1,2);
% plot(f_resampled, abs(modulated4_freq_shifted));
% title('Spectrum of Resampled Signal 4 (After resampling)');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;
% 
% 
% % --- Plot Signal 5 AFTER resampling: Time-Domain and Frequency-Domain ---
% figure;
% subplot(2,1,1);
% plot(t_new, mono5_resampled);
% title('Time-Domain Plot of Resampled Signal 5 (After resampling)');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% modulated5_freq = fft(mono5_resampled);
% modulated5_freq_shifted = fftshift(modulated5_freq);
% 
% subplot(2,1,2);
% plot(f_resampled, abs(modulated5_freq_shifted));
% title('Spectrum of Resampled Signal 5 (After resampling)');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;

BW_arr = [10e3 10e3 10e3]; % Baseband BW of signals

% AM

% Carrier frequencies for DSB-SC modulation
fc1 = 100e3; % Carrier frequency for signal 1 (100 kHz)
fc2 = 150e3; % Carrier frequency for signal 2 (100 + Δf with Δf = 50 kHz)
fc3 = 200e3; % Carrier frequency for signal 3 (100 + 2*Δf)
% fc4 = 250e3; % Carrier frequency for signal 4 (100 + 3*Δf)
% fc5 = 300e3; % Carrier frequency for signal 5 (100 + 4*Δf)


% Generate carriers
carrier1 = cos(2 * pi * fc1 * t_new); % Carrier for signal 1
carrier2 = cos(2 * pi * fc2 * t_new); % Carrier for signal 2
carrier3 = cos(2 * pi * fc3 * t_new); % Carrier for signal 3
% carrier4 = cos(2 * pi * fc4 * t_new); % Carrier for signal 4
% carrier5 = cos(2 * pi * fc5 * t_new); % Carrier for signal 5


% Modulate the signals using DSB-SC
modulated1 = mono1_resampled .* carrier1'; % Signal 1 modulated at 100 kHz
modulated2 = mono2_resampled .* carrier2'; % Signal 2 modulated at 150 kHz
modulated3 = mono3_resampled .* carrier3'; % Signal 3 modulated at 200 kHz
% modulated4 = mono4_resampled .* carrier4'; % Signal 4 modulated at 250 kHz
% modulated5 = mono5_resampled .* carrier5'; % Signal 5 modulated at 300 kHz

% Construct the FDM signal by summing the modulated signals
fdm_signal = modulated1 + modulated2 + modulated3;

% % --- Plot modulated signals: Time-Domain and Frequency-Domain ---
% figure;
% subplot(2,1,1);
% plot(t_new, modulated1);
% title('Time-Domain Plot of modulated signal 1');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % FFT for resampled mono2
% modulated2_freq = fft(modulated1);
% modulated_freq_shifted = fftshift(modulated2_freq);
% 
% subplot(2,1,2);
% plot(f_resampled, abs(modulated_freq_shifted));
% title('Spectrum of modulated signal 1 ');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;
% 
% 
% 
% figure;
% subplot(2,1,1);
% plot(t_new, modulated2);
% title('Time-Domain Plot of modulated signal 2');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % FFT for resampled mono2
% modulated2_freq = fft(modulated2);
% modulated2_freq_shifted = fftshift(modulated2_freq);
% 
% subplot(2,1,2);
% plot(f_resampled, abs(modulated2_freq_shifted));
% title('Spectrum of modulated signal 2');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;

% % --- Plot FDM signal: Time-Domain and Frequency-Domain ---
figure;
subplot(2,1,1);
plot(t_new, fdm_signal);
title('Time-Domain Plot of modulated signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT for resampled mono2
modulated_freq = fft(fdm_signal);
modulated_freq_shifted = fftshift(modulated_freq);

subplot(2,1,2);
plot(f_resampled, abs(modulated_freq_shifted));
title('Spectrum of modulated signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% RF BPF
signal_num = input("enter signal number (0 for first signal, 1 for second signal ..etc");
BandPassFilt_RF = RF_BPF(signal_num, BW_arr,new_fs);

% Apply the Bandpass Filter to FDM Signal
RF_BPF_out = filter(BandPassFilt_RF, fdm_signal);

% --- Plot Filtered Signal: Time-Domain and Frequency-Domain ---
figure;

% Time-Domain Plot
subplot(2,1,1);
plot(t_new, RF_BPF_out);
title('Time-Domain Plot of output of RF BPF ');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT of Filtered Signal
N_filtered = length(RF_BPF_out);
filtered_fft = fft(RF_BPF_out);
filtered_fft_shifted = fftshift(filtered_fft);
f_filtered = (-N_filtered/2:N_filtered/2-1) * (new_fs / N_filtered);

% Frequency-Domain Plot
subplot(2,1,2);
plot(f_filtered, abs(filtered_fft_shifted));
title('Spectrum of output of RF BPF');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% Define Intermediate Frequency
f_IF = 25e3; % Intermediate frequency (25 kHz)

% Define Oscillator Frequency
f_osc = (100e3 + signal_num * 50e3) + f_IF;

% Generate Oscillator Signal
oscillator_signal = cos(2 * pi * f_osc * t_new);

% --- Mix Filtered Signal with Oscillator Signal ---
IF_signal = RF_BPF_out .* oscillator_signal';

% % --- Plot Intermediate Frequency (IF) Signal ---
% figure;
% 
% % Time-Domain Plot
% subplot(2,1,1);
% plot(t_new, IF_signal);
% title('Time-Domain Plot of Intermediate Frequency (IF) Signal');
% xlabel('Time (s)');
% ylabel('Amplitude');
% grid on;
% 
% % FFT of IF Signal
% IF_fft = fft(IF_signal);
% IF_fft_shifted = fftshift(IF_fft);
% 
% % Frequency-Domain Plot
% subplot(2,1,2);
% plot(f_resampled, abs(IF_fft_shifted));
% title('Spectrum of Intermediate Frequency (IF) Signal');
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% grid on;



% IF BPF

BandPassFilt_IF = IF_BPF(0, BW_arr, new_fs);

% Apply the Bandpass Filter to FDM Signal
IF_BPF_out = filter(BandPassFilt_IF, IF_signal);

% --- Plot Intermediate Frequency (IF) Signal ---
figure;

% Time-Domain Plot
subplot(2,1,1);
plot(t_new, IF_BPF_out);
title('Time-Domain Plot of IF BPF output');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT of IF Signal
IF_BPF_out_fft = fft(IF_BPF_out);
IF_BPF_out_fft_shift = fftshift(IF_BPF_out_fft);

% Frequency-Domain Plot
subplot(2,1,2);
plot(f_resampled, abs(IF_BPF_out_fft_shift));
title('Spectrum of IF BPF output');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Baseband detection

% --- Mix with Local Carrier (for Baseband Detection) ---
local_carrier = cos(2 * pi * f_IF * t_new); % Local oscillator signal


baseband_mixed = IF_BPF_out .* local_carrier'; % Mix IF signal with local carrier

figure;
% Time-Domain Plot
subplot(2, 1, 1);
plot(t_new, baseband_mixed);
title('Time-Domain Plot of Baseband mixed Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT of Baseband Signal
baseband_mixed_fft = fft(baseband_mixed);
baseband_mixed_fft_shift = fftshift(baseband_mixed_fft);

% Frequency-Domain Plot
subplot(2, 1, 2);
plot(f_resampled, abs(baseband_mixed_fft_shift));
title('Spectrum of Baseband mixed Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% --- Define LPF Specifications ---
F_stop = 20e3;  % Stopband frequency = 20 kHz
F_pass = 10e3;  % Passband frequency = 10 kHz
A_stop = 60;    % Attenuation in the stopband = 60 dB
A_pass = 1;     % Amount of ripple allowed in the passband = 1 dB

% Create Low-Pass Filter Specification Object
LPFSpecObj = fdesign.lowpass('Fp,Fst,Ap,Ast', F_pass, F_stop, A_pass, A_stop, new_fs);

% Design the LPF using Butterworth Filter
LPF = design(LPFSpecObj, 'butter');

% --- Visualize the Low-Pass Filter Characteristics ---
%fvtool(LPF); % to view LPF response


% --- Apply Low-Pass Filter to Extract Baseband Signal ---
baseband_signal = filter(LPF, baseband_mixed);

% --- Plot Baseband Signal: Time-Domain and Frequency-Domain ---
figure;

% Time-Domain Plot
subplot(2, 1, 1);
plot(t_new, baseband_signal);
title('Time-Domain Plot of Baseband Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% FFT of Baseband Signal
baseband_fft = fft(baseband_signal);
baseband_fft_shifted = fftshift(baseband_fft);

% Frequency-Domain Plot
subplot(2, 1, 2);
plot(f_resampled, abs(baseband_fft_shifted));
title('Spectrum of Baseband Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;


% Normalize the baseband signal
baseband_signal_normalized = baseband_signal / max(abs(baseband_signal));


% Downsample to 44.1 kHz
fs_downsampled = 44100; % Target sampling frequency
downsample_factor = new_fs / fs_downsampled;
baseband_downsampled = downsample(baseband_signal_normalized, round(downsample_factor));

% Play the downsampled signal
sound(baseband_downsampled, fs_downsampled);

