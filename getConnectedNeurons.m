function connectedNeuronIndices = getConnectedNeurons(network)
    connectedNeuronIndices = [];
    for i = 1:length(network)
        if any(network(i, :))
            connectedNeuronIndices = [connectedNeuronIndices i];
        end
    end
end