function corrCoeffMatrix = getBinSpikesCorrCoeffs(binnedSpikeCell)
    corrCoeffMatrix = NaN(length(binnedSpikeCell), length(binnedSpikeCell));
    for neuron = 1:length(binnedSpikeCell)
        for compNeuron = 1:length(binnedSpikeCell)
            oneNeuron = binnedSpikeCell{neuron};
            otherNeuron = binnedSpikeCell{compNeuron};
            if isempty(oneNeuron) || isempty(otherNeuron)
                continue
            end
                corrVal = corrcoef(double(binnedSpikeCell{neuron}), double(otherNeuron(1:length(oneNeuron))));
                corrCoeffMatrix(neuron, compNeuron) = corrVal(1,2);
        end    
    end
    %corrCoeffMatrix = corrCoeffMatrix(~isnan(corrCoeffMatrix));
end

