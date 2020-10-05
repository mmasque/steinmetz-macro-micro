function spikeTimesCell = reformatSpikeTimes(indices, times, clusters)
 
    spikeTimesCell = cell(1, max(indices));
    sortedClustersAndTimes = sortrows([double(clusters) times]);
    %naive: for i in indices, get times at the appropriate clusters. 
    for ind = indices'
       indices_of_i = ismember(sortedClustersAndTimes(:, 1), ind);
       spikeTimesCell{ind} = (sortedClustersAndTimes(indices_of_i, 2));
    end
    
end

    