function simSpikesMat = allSimSpikes(binnedSpikeTimesCell)
    simSpikesMat = NaN(length(binnedSpikeTimesCell),length(binnedSpikeTimesCell));
    for i = 1:length(binnedSpikeTimesCell)
        if isempty(binnedSpikeTimesCell{i})
            continue
        end
        for j = 1:length(binnedSpikeTimesCell)
            if isempty(binnedSpikeTimesCell{j})
                continue
            end
            simSpikesMat(i,j) = countSimultaneousSpikes(binnedSpikeTimesCell{i}, binnedSpikeTimesCell{j});
        end
    end
end