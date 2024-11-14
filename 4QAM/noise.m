function [narrowband_noises] = noise(sigma, t)

    % Initial setup for generating narrowband noise
    Fs = 100e3;  % Sampling frequency (100 kHz)

    % Prepare an array to store each narrowband noise signal
    narrowband_noises = zeros(length(sigma), length(t));

    % Generate narrowband noise for each value of sigma (noise variance)
    for i = 1:length(sigma)
        % Create Gaussian random noise with mean 0 and standard deviation based on sigma(i)
        gaussian_noise = sigma(i) * randn(size(t));

        % Set the carrier frequency and bandwidth for the noise
        f_center = 20e3;  % Center frequency of the noise (20 kHz)
        bandwidth = 3.2e3;  % Bandwidth of the noise (3.2 kHz)

        % Design a bandpass filter to create narrowband noise
        d = designfilt('bandpassfir', 'FilterOrder', 50, ...  % Filter order of 50
            'CutoffFrequency1', f_center-bandwidth/2, ...    % Lower cutoff frequency
            'CutoffFrequency2', f_center+bandwidth/2, ...    % Upper cutoff frequency
            'SampleRate', Fs);  % Sampling frequency

        % Filter the Gaussian noise to create narrowband noise
        narrowband_noises(i, :) = filter(d, gaussian_noise);
    end
end
