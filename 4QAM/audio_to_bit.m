function second_database = audio_to_bit(audio_file)
    % Load the audio file
    [record, fs_record] = audioread(audio_file);

    % Normalize the signal to the range [0, 1]
    min_record = min(record);
    max_record = max(record);
    normalized_record = (record - min_record) / (max_record - min_record);

    % Calculate the downsampling factor corresponding to 1.25 milliseconds
    sampling_interval = 1.25e-3; % 1.25 milliseconds
    samples_per_interval = round(fs_record * sampling_interval);

    % Perform downsampling every 1.25 milliseconds
    sampled_record = downsample(normalized_record, samples_per_interval);

    % Perform quantization to 16 levels (4 bits)
    NQ = 16; % Number of quantization levels
    quantized_record = round(sampled_record * (NQ - 1));

    % Check the number of quantization levels used
    actual_levels = numel(unique(quantized_record));
    disp(['Number of quantization levels used: ', num2str(actual_levels)]);

    % Convert the quantized values to binary representation
    bit_depth = 4; % Number of bits
    binary_record = dec2bin(quantized_record, bit_depth) - '0'; % Direct conversion to numerical array

    % Reshape to 2 columns
    second_database = reshape(binary_record', 2, [])';

    % Calculate the new sampling frequency after downsampling
    new_fs = fs_record / samples_per_interval;
    disp(['New sampling frequency: ', num2str(new_fs), ' Hz']);
end
