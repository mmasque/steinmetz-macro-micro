function binarySpikes = spikeTimesToBinary(binnedSpikeTimesCell)
    binarySpikes = cell(size(binnedSpikeTimesCell));
    for c = 1:length(binnedSpikeTimesCell)
        curr = binnedSpikeTimesCell{c};
        if isempty(curr)
            continue
        end
        binarySpikes{c} = curr > 0;
    end
end