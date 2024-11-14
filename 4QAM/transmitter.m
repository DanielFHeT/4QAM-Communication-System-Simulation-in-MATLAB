function s_m = transmitter(a_n, b_n, omega_c, a_c, t)
% Repeat the symbols efficiently using repelem
% a_n_repeated and b_n_repeated extend the real and imaginary symbol sequences to match the length of time vector t
a_n_repeated = repelem(a_n, length(t) / length(a_n));
b_n_repeated = repelem(b_n, length(t) / length(b_n));

% Create the cosine and sine carriers
I_carrier = cos(omega_c * t); % In-phase (I) carrier: cosine wave with frequency omega_c
Q_carrier = sin(omega_c * t); % Quadrature (Q) carrier: sine wave with frequency omega_c

% Modulation
% Perform IQ modulation where the signal is a combination of the real and imaginary parts
% a_n_repeated is modulated with the I_carrier and b_n_repeated is modulated with the Q_carrier
s_m = a_c * (a_n_repeated .* I_carrier - b_n_repeated .* Q_carrier); % Final modulated signal
end
