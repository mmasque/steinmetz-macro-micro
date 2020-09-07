function [indices, clusterIDs, times] = getGoodIndices(path)

    phy_annotation = readNPY(strcat(path, "/clusters._phy_annotation.npy"));
    indices = find(phy_annotation >= 2);


    spike_clusters = readNPY(strcat(path, "/spikes.clusters.npy"));
    spike_times = readNPY(strcat(path, "/spikes.times.npy"));
    
    all_indices = find(ismember(spike_clusters, indices));
    clusterIDs = spike_clusters(all_indices);
    times = spike_times(all_indices);
    
end
