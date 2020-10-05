function connectedNeuronIndices = getConnectedNeurons(network)
    connectedNeuronIndices = [];
    for i = 1:length(network)
        if any(network(i, :))
            connections = find(network(i,:))';
            indexpairs = [ones(length(connections), 1) * i, connections];
            connectedNeuronIndices = [connectedNeuronIndices;indexpairs];%[connectedNeuronIndices i];
        end
    end
end