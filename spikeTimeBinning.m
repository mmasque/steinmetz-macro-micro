function allCounts = spikeTimeBinning(spikeTimesCell, binWidth)
    %{
        Counts spikes in each bin
        - binWidth: in seconds 
    %}
    [~, columns] = size(spikeTimesCell); % Get size of cell array
    allMin = Inf;
    allMax = -Inf;
    allCounts = cell(size(spikeTimesCell));
    for col = 1 : columns
        theseTimes = spikeTimesCell{1, col};
        if isempty(theseTimes)
            continue;
        end
        thisMin = min(theseTimes);
        thisMax = max(theseTimes);
        
        allMin = min(allMin, thisMin);
        allMax = max(allMax, thisMax);
    end 
    % Now construct edges
    edges = allMin : binWidth : (allMax + binWidth);
    
    for col = 1 : columns
        theseTimes = spikeTimesCell{1, col};
        allCounts{col} = (histcounts(theseTimes, edges));
    end
     
end