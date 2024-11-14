function plotConstellation(a_n, b_n, A_n, constellation_name, letter, scale, isTransmitter)
    % Open a new figure for each constellation plot
    figure;
    
    % Plot the constellation points
    scatter(a_n, b_n, 'filled');
    xlabel('In-phase Component');
    ylabel('Quadrature Component');
    title(constellation_name);
    grid on;
    axis([-3 3 -3 3] * scale);
    
    hold on; % Allow annotations or further plotting
    
    if isTransmitter
        % Define the Gray code mapping (as a string for annotation purposes)
        gray_code_bits = ["00", "01", "10", "11"];
        
        % Annotate the constellation points with Gray code and labels
        for i = 1:length(a_n)
            text(a_n(i) + 0.05 * scale, b_n(i) + 0.05 * scale, ...
                sprintf('A%d (%s) [%s]', i, string(A_n(i)), gray_code_bits(i)), 'FontSize', 8);
        end
    else
        % Label points on the plot with simple labels (without Gray code)
        for i = 1:length(a_n)
            text(a_n(i) + 0.05 * scale, b_n(i), [letter, num2str(i)], 'FontSize', 12);
        end

        % Compute and display distances between each symbol
        fprintf('Distances between symbols in %s:\n', constellation_name);
        for i = 1:length(a_n)
            for j = i+1:length(a_n)
                distance = sqrt((a_n(i) - a_n(j))^2 + (b_n(i) - b_n(j))^2);
                fprintf('Distance between %s%d and %s%d: %.5f\n', letter, i, letter, j, distance);
            end
        end
    end
    
    hold off; % Finish annotation or plotting
end