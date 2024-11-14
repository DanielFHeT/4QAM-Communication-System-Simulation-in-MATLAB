function [qki, qkq] = receiver(s_r, a_0, pshi, omega_c, samples, len_a_n, t)
    % Function to calculate qki and qkq based on the input signal and parameters

    % Reshape the time vector t and the received signal s_r into matrices with rows corresponding to symbols
    t_matrix = reshape(t, samples, len_a_n)'; % Reshape the time vector to match the number of symbols (len_a_n)
    s_m_matrix = reshape(s_r, samples, len_a_n)'; % Reshape the received signal matrix based on symbol length

    % Calculate the integrand parameters for qki and qkq

    % Compute the integrand for qki using the cosine term
    integrand_qki = s_m_matrix .* a_0 .* cos(omega_c .* t_matrix + pshi); % For qki, phase shift (pshi) included

    % Compute the integrand for qkq using the sine term with a negative sign
    integrand_qkq = -s_m_matrix .* a_0 .* sin(omega_c .* t_matrix + pshi); % For qkq, phase shift (pshi) included

    % Calculate the time step differences (dt) between consecutive time points
    dt = diff(t_matrix, 1, 2);

    % Perform trapezoidal integration manually for qki and qkq over time
    % The integration is done symbol by symbol using the trapezoidal rule
    qki = sum((integrand_qki(:, 1:end-1) + integrand_qki(:, 2:end)) .* dt / 2, 2); % Trapezoidal integration for qki
    qkq = sum((integrand_qkq(:, 1:end-1) + integrand_qkq(:, 2:end)) .* dt / 2, 2); % Trapezoidal integration for qkq

    % Transpose qki and qkq to return them as row vectors
    qki = qki'; 
    qkq = qkq';
end
