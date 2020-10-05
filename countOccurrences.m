function occurrences = countOccurrences(data, nValues, nStates)
    % data should be nmeasurements x nValues (number of 'nodes')
    nPoss = nStates ^ nValues;
    occurrences = zeros(1, nPoss);

    for i = 1:nPoss
        state_sample = loli_index2state(i, nStates);
        state_sample = logical(state_sample);
        state_match = ismember(data, state_sample, 'rows');
        state_positions = find(state_match);
        
        occurrences(i) = length(state_positions);
    end
end
