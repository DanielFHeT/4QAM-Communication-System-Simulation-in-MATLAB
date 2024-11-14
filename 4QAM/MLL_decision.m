function [close_symbol, distances] = MLL_decision(received_symbols, constellation)
    % Calculate the distances between each received symbol and every point in the constellation
    % The abs function computes the Euclidean distance (since symbols are complex numbers)
    distances = abs(received_symbols.' - constellation);
    
    % Find the closest constellation point for each received symbol
    % min function returns the index of the closest symbol in the constellation
    [~, closest_symbol_indices] = min(distances, [], 2);
    
    % Return the MLL decision - the closest point in the constellation for each received symbol
    close_symbol = constellation(closest_symbol_indices);
end
