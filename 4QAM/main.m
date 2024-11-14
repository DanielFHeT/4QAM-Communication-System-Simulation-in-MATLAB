clc;
clear;
% Simulation Chapter 3

%% step 1 crate database

% Section 1a - Gray coding map
First_database = [0, 0; 0, 1; 1, 0; 1, 1];

% Section 1b - Generate random bit array for the base dataset
second_database = audio_to_bit('record.opus');

%% step 2 convert bits to symbol
[A_n1, a_n1, b_n1] = Symbol_Encoder(First_database);

disp('Encoded symbols from First_database:');
disp(A_n1);

plotConstellation(a_n1, b_n1, A_n1, 'Transmitter Constellation', 'A', 1, true);


% convert bits to symbol the second_database
[A_n2, a_n2, b_n2] = Symbol_Encoder(second_database);

disp('First 10 encoded symbols from second_database:');
disp(A_n2(1:10));  % Display the first 10 encoded symbols for brevity

%% step 3 transmitter
% Modulation on First_database symbols

P_c = 10; % Carrier power

a_c = sqrt(2 * P_c); % Amplitude of the carrier signal

T_s = 0.625e-3; % Symbol period (time for one symbol)

f_c = 20e3; % Carrier frequency (Hertz)

omega_c = 2 * pi * f_c; % Angular frequency of the carrier

samples = 200; % Number of samples per symbol

% Generate a time vector for sampling, with enough points to cover all symbols
t1 = linspace(0, length(A_n1)*T_s, length(A_n1)*samples); 

% Modulate the input symbol sequences a_n1 and b_n1 using the transmitter function
s_m1 = transmitter(a_n1,b_n1,omega_c,a_c,t1);

% Plot the modulated signal
figure;
plot(t1, s_m1);
xlabel("Time (s)"); % Label the x-axis as time in seconds
ylabel("Amplitude"); % Label the y-axis as signal amplitude
title("Modulated Signal for First Database") % Set the title for the plot
grid on; % Turn on the grid for better readability







%% 
% Modulation on second database symbols followed by frequency conversion

t2 = linspace(0, length(A_n2)*T_s, length(A_n2)*samples); % Time vector for sampling, enough points to cover all symbols

s_m2 = transmitter(a_n2,b_n2,omega_c,a_c,t2); % Modulate the input symbol sequences a_n2 and b_n2 using the transmitter function

Y = fftshift(fft(s_m2)); % Compute the Fourier Transform of the modulated signal and shift it to center zero frequency

% Normalize the FFT to get the magnitude of the frequency components
Y_normalized = abs(Y) / length(Y);

% Create a frequency vector to cover both positive and negative frequencies
f = linspace(-0.5, 0.5, length(Y)) * (1 / T_s) * 200; % Frequency vector scaled to match the sampling rate

% Plot the normalized frequency spectrum
figure;
plot(f, Y_normalized); % Plot the frequency magnitude against the frequency vector
title('Frequency Spectrum of the Modulated Signal'); % Set title for the plot
xlabel('Frequency (Hz)'); % Label the x-axis as frequency in Hz
ylabel('Magnitude'); % Label the y-axis as magnitude of the frequency components
grid on; % Turn on the grid for better readability
xlim([-9.8e4 9.8e4]); % Set the x-axis limits for the frequency range

%% receiver & MLL decision & symbol_decoder with pshi = 0

% receiver 
a_0 = 1; % Amplitude of the received signal

pshi1 = 0; % Initial phase shift of the received signal

len_A_n = length(A_n1); % Length of the first database (number of transmitted symbols)

% Receiver function to demodulate the received signal s_m1
% qki and qkq are the in-phase (I) and quadrature (Q) components of the received signal
[qki, qkq] = receiver(s_m1, a_0, pshi1, omega_c, samples, len_A_n, t1);

% Combine the in-phase and quadrature components to form the complex received symbols
received_symbols_1 = qki + 1j*qkq;

% MLL decision 
% Compute the closest transmitted symbols (MLL decision) based on the received symbols
% close_symbol represents the estimated transmitted symbols
% distances represent the Euclidean distances between the received symbols and each possible constellation point
[close_symbol, distances] = MLL_decision(received_symbols_1, A_n1);

% Display the results
disp('pshi = 0:')

disp('MLL_decision:')
disp(close_symbol) % Display the closest symbol decisions

disp('distances:')
disp(distances) % Display the distances from the received symbols to the closest symbols

% symbol_decoder
% Decode the symbols back into the original bit sequence using the symbol_decoder function
b_r = symbol_decoder(close_symbol);
disp('symbol_decoder:');
disp(b_r) % Display the decoded bit sequence



%% receiver on First_database symbols with pshi = 30
% receiver 
pshi2 = deg2rad(30); % Convert 30 degrees to radians for the phase shift

% Receiver function to demodulate the received signal s_m1 with phase shift pshi2
% qki and qkq are the in-phase (I) and quadrature (Q) components of the received signal
[qki, qkq] = receiver(s_m1, a_0, pshi2, omega_c, samples, len_A_n, t1);

% Combine the in-phase and quadrature components to form the complex received symbols
received_symbols_2 = qki + 1j*qkq;

% MLL decision 
% Compute the closest transmitted symbols (MLL decision) based on the received symbols
% close_symbol represents the estimated transmitted symbols with the phase shift applied
% distances represent the Euclidean distances between the received symbols and each possible constellation point
[close_symbol, distances] = MLL_decision(received_symbols_2, A_n1);

% Display the results for phase shift of 30 degrees
disp('pshi = 30:')

disp('MLL_decision:')
disp(close_symbol) % Display the closest symbol decisions

disp('distances:')
disp(distances) % Display the distances from the received symbols to the closest symbols

% symbol_decoder
% Decode the symbols back into the original bit sequence using the symbol_decoder function
b_r = symbol_decoder(close_symbol);
disp('symbol_decoder:');
disp(b_r) % Display the decoded bit sequence


%%
% signal Before multiplication by cos(2 * pi * f_c * t) sin(2 * pi * f_c * t) in the transmitter on First_database symbols

% Repeat the real (a_n1) and imaginary (b_n1) parts of the symbols over the number of samples for plotting
a_n_repeated = repelem(a_n1, samples);
b_n_repeated = repelem(b_n1, samples);

% Plot the signal before multiplying by cos(2 * pi * f_c * t) and sin(2 * pi * f_c * t) in the transmitter
figure;
plot(t1, a_n_repeated, 'b-', 'LineWidth', 2); % Blue solid line for real part of the signal
hold on;
plot(t1, b_n_repeated, 'r-', 'LineWidth', 2, 'Color', [1 0 0 0.5]); % Red transparent line for imaginary part of the signal
title('Modulated Signal for First Database before multiplying by sine cosine');
xlabel('Time (s)'); % Label for x-axis: Time in seconds
ylabel('Amplitude'); % Label for y-axis: Amplitude of the signal
legend('Real','Imaginary','Fontsize',10) % Legend for real and imaginary parts
grid on; % Enable grid for better readability
hold off;

% signal after match filter in receiver with pshi=0 on First_database symbols

% Repeat the real and imaginary parts of the received symbols after matched filtering (pshi = 0)
real_pshi_0 = repelem(real(received_symbols_1), samples);
image_pshi_0 = repelem(imag(received_symbols_1), samples);

% Plot the received signal after matched filtering with pshi = 0
figure();
plot(t1, real_pshi_0, 'b-', 'LineWidth', 2); % Blue solid line for real part of the received signal
hold on;
plot(t1, image_pshi_0, 'r-', 'LineWidth', 2, 'Color', [1 0 0 0.5]); % Red transparent line for imaginary part of the received signal
title('Received Signal for First Dataset After Matched Filtering with pshi=0');
xlabel('Time (s)'); % Label for x-axis: Time in seconds
ylabel('Amplitude'); % Label for y-axis: Amplitude of the received signal
legend('Real','Imaginary','Fontsize',10) % Legend for real and imaginary parts
grid on; % Enable grid for better readability
hold off;

% signal after match filter in receiver with pshi=30 on First_database symbols

% Repeat the real and imaginary parts of the received symbols after matched filtering (pshi = 30)
real_pshi_30 = repelem(real(received_symbols_2), samples);
image_pshi_30 = repelem(imag(received_symbols_2), samples);

% Plot the received signal after matched filtering with pshi = 30
figure();
plot(t1, real_pshi_30, 'b-', 'LineWidth', 2); % Blue solid line for real part of the received signal
hold on;
plot(t1, image_pshi_30, 'r-', 'LineWidth', 2, 'Color', [1 0 0 0.5]); % Red transparent line for imaginary part of the received signal
title('Received Signal for First Dataset After Matched Filtering with pshi=30');
xlabel('Time (s)'); % Label for x-axis: Time in seconds
ylabel('Amplitude'); % Label for y-axis: Amplitude of the received signal
legend('Real','Imaginary','Fontsize',10) % Legend for real and imaginary parts
grid on; % Enable grid for better readability
hold off;



%% constellation and Distances between symbols

% Transmitter constellation
% Plot the transmitter constellation using the a_n1 and b_n1 (real and imaginary parts) and A_n1 (symbol annotations)
% The title of the plot will be 'Transmitter Constellation'
% 'A' is the prefix for annotations, 1 is the scaling factor for symbol distances, and false means no extra formatting
plotConstellation(a_n1, b_n1, A_n1, 'Transmitter Constellation', 'A', 1, false);

% Receiver Constellation with pshi = 0
% Plot the receiver constellation for received_symbols_1 (pshi = 0)
% Only the real and imaginary parts of the received symbols are plotted without annotations
% The title of the plot will be 'Receiver Constellation with pshi = 0'
% 'a' is the prefix for annotations, 1e-3 is the scaling factor for symbol distances, and false means no extra formatting
plotConstellation(real(received_symbols_1), imag(received_symbols_1), [], 'Receiver Constellation with pshi = 0', 'a', 1e-3, false);

% Receiver Constellation with pshi = 30
% Plot the receiver constellation for received_symbols_2 (pshi = 30)
% Only the real and imaginary parts of the received symbols are plotted without annotations
% The title of the plot will be 'Receiver Constellation with pshi = 30'
% 'a' is the prefix for annotations, 1e-3 is the scaling factor for symbol distances, and false means no extra formatting
plotConstellation(real(received_symbols_2), imag(received_symbols_2), [], 'Receiver Constellation with pshi = 30', 'a', 1e-3, false);

%% chapter 4
% Parameter definitions
A_0 = 1.32e6; % Amplitude of the transmitted signal
P_r = P_c;  % Received power, assumed equal to transmitted power (10W)
E_g = T_s;  % Symbol duration (0.625 ms)

% Decimal SNR values
per_symbol_max = 1e-3; % Maximum symbol error rate
per_symbol_min = 2e-1; % Minimum symbol error rate

% Decimal calculation of gamma_d for the two cases
gamma_d_max = 2 * (erfcinv(per_symbol_max))^2; % Maximum gamma_d based on error rate
gamma_d_min = 2 * (erfcinv(per_symbol_min))^2; % Minimum gamma_d based on error rate

Gamma_d_symbol = gamma_d_min:2:gamma_d_max; % SNR (Gamma_d) values for each symbol

% Compute N_0 (noise power spectral density) for each Gamma_d_symbol value
N_0 = (2 * P_r * E_g) ./ Gamma_d_symbol;

% Calculate the noise variance (sigma) for each value of N_0
sigma = sqrt((A_0^2 .* N_0 .* E_g) / 4);

% Generate narrowband noise for the signal
narrowband_noises = noise(sigma, t2);

% Add noise to the signal
x_t = s_m2 + narrowband_noises;

% Get dimensions of x_t
[m, n] = size(x_t); % m is the number of signals, n is the number of samples

% Length of the second database
len_A_n2 = length(A_n2);

% Define variables for qki and qkq (real and imaginary parts of received symbols)
qki = zeros(m, len_A_n2); % Create matrix for qki (real part)
qkq = zeros(m, len_A_n2); % Create matrix for qkq (imaginary part)

% Assumption: pshi=0 (no phase shift)

% Calculate qki and qkq for each row of x_t (for each noisy signal)
for i = 1:m
    [qki(i, :), qkq(i, :)] = receiver(x_t(i, :), A_0, pshi1, omega_c, samples, len_A_n2, t2);
end

% Compute received_symbols_3 (complex received symbols)
received_symbols_3 = qki + 1j * qkq;

% Calculate errors and BER (Bit Error Rate)
errors = zeros(1, length(m)); % Initialize error vector
ber = zeros(1, length(m)); % Initialize BER vector
num_bits_per_symbol = 2; % For 4QAM, there are 2 bits per symbol

for i = 1:m
    for j = 1:len_A_n2
        [close_symbol(i, j), min_distances] = MLL_decision(received_symbols_3(i, j), A_n1); % Perform MLL decision
    end
    errors(i) = sum(close_symbol(i, :) ~= A_n2);  % Total errors for each symbol
    ber(i) = errors(i) / len_A_n2;  % Symbol error rate (BER for symbols)
end

total_bits = len_A_n2 * num_bits_per_symbol;  % Total number of bits transmitted
ber_bit = errors / total_bits;  % Bit Error Rate (BER) for each iteration

ber_bit_Full_sync = ber_bit; % Save BER values for full synchronization

M = 4; % Modulation order (4 for 4QAM)

% Theoretical calculation of BER
PerSym = erfc(sqrt(Gamma_d_symbol) / sqrt(2)); % Theoretical symbol error rate
PerBit = PerSym / log2(M); % Theoretical bit error rate

% Convert Gamma_d_symbol to dB for symbols
Gamma_d_symbol_dB = 10 * log10(Gamma_d_symbol);

% Convert Gamma_d_symbol to dB for bits
Gamma_d_bit_dB = Gamma_d_symbol_dB / log2(M);

% Plot simulation results only
figure("Name", 'Simulation Only', 'NumberTitle', 'off');
semilogy(Gamma_d_bit_dB, ber_bit, 'LineWidth', 2.5, 'DisplayName', 'Simulation');
xlabel('SNR [dB]'); % X-axis label: Signal-to-Noise Ratio in dB
ylabel('Bit Error Probability'); % Y-axis label: Bit Error Probability
grid on; % Enable grid for better readability
legend('show'); % Display legend
title('Simulation BER vs SNR'); % Title of the plot
hold off;

% Plot comparison between simulation and theoretical BER
figure("Name", 'Simulation Vs Theoretical', 'NumberTitle', 'off');
semilogy(Gamma_d_bit_dB, ber_bit, 'LineWidth', 2.5, 'DisplayName', 'Simulation'); % Plot simulation results
hold on;
semilogy(Gamma_d_bit_dB, PerBit, 'LineWidth', 2.5, 'DisplayName', 'Theoretical'); % Plot theoretical results
xlabel('SNR [dB]'); % X-axis label: Signal-to-Noise Ratio in dB
ylabel('Bit Error Probability'); % Y-axis label: Bit Error Probability
grid on; % Enable grid for better readability
legend('show'); % Display legend
title('Simulation vs Theoretical BER'); % Title of the plot
hold off;


%% Chapter 5 - Power Penalty Calculation for Different Phase Deviations

% Calculate power penalty for a phase deviation of 5 degrees
qki = zeros(m, len_A_n2); % Create matrix for qki (real part)
qkq = zeros(m, len_A_n2); % Create matrix for qkq (imaginary part)

pshi2 = deg2rad(5); % Convert 5 degrees to radians for phase deviation

% Calculate qki and qkq for each row of x_t (for each noisy signal)
for i = 1:m
    [qki(i, :), qkq(i, :)] = receiver(x_t(i, :), A_0, pshi2, omega_c, samples, len_A_n2, t2);
end

% Compute received_symbols_3 (complex received symbols)
received_symbols_3 = qki + 1j * qkq;

% Calculate errors and BER
errors = zeros(1, m); % Initialize errors vector
ber = zeros(1, m); % Initialize BER vector
num_bits_per_symbol = 2; % For 4QAM, there are 2 bits per symbol

for i = 1:m
    for j = 1:len_A_n2
        [close_symbol(i, j), min_distances] = MLL_decision(received_symbols_3(i, j), A_n1); % Perform MLL decision
    end
    errors(i) = sum(close_symbol(i, :) ~= A_n2);  % Total errors for each symbol
    ber(i) = errors(i) / len_A_n2;  % Symbol error rate (BER for symbols)
end

total_bits = len_A_n2 * num_bits_per_symbol;  % Total number of bits transmitted
ber_bit_phase_deviation_5 = errors / total_bits;  % Bit Error Rate (BER) for each iteration with 5-degree phase deviation

% Plot Simulation vs Theoretical for psi = 5 degrees
figure("Name", 'Simulation Vs Theoretical psi = 5 degrees', 'NumberTitle', 'off');
semilogy(Gamma_d_bit_dB, ber_bit_phase_deviation_5, 'LineWidth', 2.5, 'DisplayName', 'Simulation psi 5');
hold on;
semilogy(Gamma_d_bit_dB, PerBit, 'LineWidth', 2.5, 'DisplayName', 'Theoretical'); % Plot theoretical BER
xlabel('SNR [dB]'); % X-axis label: Signal-to-Noise Ratio in dB
ylabel('Bit Error Probability'); % Y-axis label: Bit Error Probability
grid on; % Enable grid for better readability
legend('show'); % Show legend

% Calculate power penalty for psi = 5 degrees
y_target = 10^-2; % Target for error rate 10^-2
x1_psi5 = interp1(ber_bit_phase_deviation_5, Gamma_d_bit_dB, y_target); % SNR for 5-degree deviation
x2_psi5 = interp1(PerBit, Gamma_d_bit_dB, y_target); % SNR for theoretical BER
distance5 = abs(x1_psi5 - x2_psi5); % Calculate the power penalty in dB
plot([x1_psi5 x2_psi5], [y_target y_target], '-'); % Plot horizontal line indicating the penalty
text((x1_psi5 + x2_psi5)/2, y_target, sprintf('Distance = %.2f dB', distance5), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off;

% Calculate power penalty for a phase deviation of 10 degrees
qki = zeros(m, len_A_n2); % Create matrix for qki (real part)
qkq = zeros(m, len_A_n2); % Create matrix for qkq (imaginary part)

pshi2 = deg2rad(10); % Convert 10 degrees to radians for phase deviation

% Calculate qki and qkq for each row of x_t (for each noisy signal)
for i = 1:m
    [qki(i, :), qkq(i, :)] = receiver(x_t(i, :), A_0, pshi2, omega_c, samples, len_A_n2, t2);
end

% Compute received_symbols_3 (complex received symbols)
received_symbols_3 = qki + 1j * qkq;

% Calculate errors and BER
errors = zeros(1, m); % Initialize errors vector
ber = zeros(1, m); % Initialize BER vector

for i = 1:m
    for j = 1:len_A_n2
        [close_symbol(i, j), min_distances] = MLL_decision(received_symbols_3(i, j), A_n1); % Perform MLL decision
    end
    errors(i) = sum(close_symbol(i, :) ~= A_n2);  % Total errors for each symbol
    ber(i) = errors(i) / len_A_n2;  % Symbol error rate (BER for symbols)
end

total_bits = len_A_n2 * num_bits_per_symbol;  % Total number of bits transmitted
ber_bit_phase_deviation_10 = errors / total_bits;  % Bit Error Rate (BER) for each iteration with 10-degree phase deviation

% Plot Simulation vs Theoretical for psi = 10 degrees
figure("Name", 'Simulation Vs Theoretical psi = 10 degrees', 'NumberTitle', 'off');
semilogy(Gamma_d_bit_dB, ber_bit_phase_deviation_10, 'LineWidth', 2.5, 'DisplayName', 'Simulation psi 10');
hold on;
semilogy(Gamma_d_bit_dB, PerBit, 'LineWidth', 2.5, 'DisplayName', 'Theoretical'); % Plot theoretical BER
xlabel('SNR [dB]'); % X-axis label: Signal-to-Noise Ratio in dB
ylabel('Bit Error Probability'); % Y-axis label: Bit Error Probability
grid on; % Enable grid for better readability
legend('show'); % Show legend

% Calculate power penalty for psi = 10 degrees
y_target = 10^-2; % Target for error rate 10^-2
x1_psi10 = interp1(ber_bit_phase_deviation_10, Gamma_d_bit_dB, y_target); % SNR for 10-degree deviation
x2_psi10 = interp1(PerBit, Gamma_d_bit_dB, y_target); % SNR for theoretical BER
distance10 = abs(x1_psi10 - x2_psi10); % Calculate the power penalty in dB
plot([x1_psi10 x2_psi10], [y_target y_target], '-'); % Plot horizontal line indicating the penalty
text((x1_psi10 + x2_psi10)/2, y_target, sprintf('Distance = %.2f dB', distance10), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
hold off;

% Convert the distances to linear scale
distance5_linear = 10.^(distance5/10); % Convert 5-degree power penalty to linear scale
distance10_linear = 10.^(distance10/10); % Convert 10-degree power penalty to linear scale

% Display the results in the console
fprintf('The power penalty of psi 5 degrees is: %f\n', distance5_linear); % Print power penalty for 5 degrees
fprintf('The power penalty of psi 10 degrees is: %f\n', distance10_linear); % Print power penalty for 10 degrees



