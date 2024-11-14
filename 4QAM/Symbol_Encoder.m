function [A_n, a_n, b_n] = Symbol_Encoder(bits)
    % Define the Gray code mapping (lookup table)
    Symbol_Gray = [1+1j, 1-1j, -1+1j, -1-1j];

    % Pre-calculate the decimal indices for all bit pairs
    decimal_indices = bits * [2; 1];  % Multiply bit columns directly for indexing

    % Map decimal indices to complex numbers using the lookup table
    A_n = Symbol_Gray(decimal_indices + 1);  % MATLAB indexing starts at 1

    % Extract real and imaginary parts
    a_n = real(A_n);  % Real part from the encoded symbols
    b_n = imag(A_n);  % Imaginary part from the encoded symbols
end
