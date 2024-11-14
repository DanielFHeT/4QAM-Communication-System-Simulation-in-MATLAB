function b_r = symbol_decoder(close_symbol)
    % Initialize the bits array with preallocation to improve performance
    % The bit array length is twice the number of received symbols, since 4QAM has 2 bits per symbol
    b_r = zeros(1, 2*length(close_symbol));

    % Loop through each received symbol to map it to its corresponding bits using Gray code
    for k = 1:length(close_symbol)
        % Extract the real and imaginary parts of the current symbol and round them to the nearest integer
        real_part = round(real(close_symbol(k)));
        imag_part = round(imag(close_symbol(k)));

        % Map the real and imaginary parts to the corresponding bit pair using Gray code
        if real_part == 1 && imag_part == 1
            b_r(2*k-1:2*k) = [0 0];  % 1+1i corresponds to bits 00
        elseif real_part == 1 && imag_part == -1
            b_r(2*k-1:2*k) = [0 1];  % 1-1i corresponds to bits 01
        elseif real_part == -1 && imag_part == 1
            b_r(2*k-1:2*k) = [1 0];  % -1+1i corresponds to bits 10
        elseif real_part == -1 && imag_part == -1
            b_r(2*k-1:2*k) = [1 1];  % -1-1i corresponds to bits 11
        end
    end
end
