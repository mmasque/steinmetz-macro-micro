path = "allData/Cori_2016-12-14/spikes.times.npy";
spike_times = readNPY(path);

path = "allData/Cori_2016-12-14/spikes.clusters.npy";
spike_clusters = readNPY(path);


path = "allData/Cori_2016-12-14/channels.brainLocation.tsv";
positions = readtable(path, 'FileType', 'text', 'Delimiter', '\t');
%positions = table2array(positions(:, 1:3));

path = "allData/Cori_2016-12-14/clusters.depths.npy";
depths = readNPY(path);

% filter by phy annotation: >= 2 is good
path = "allData/Cori_2016-12-14/clusters._phy_annotation.npy";
phy_annotation = readNPY(path);

good_cluster_indices = find(phy_annotation >= 2);

% find the spikes that are good
good_spikes = ismember(spike_clusters, good_cluster_indices);   % binary search! 
spike_indices = find(good_spikes);
% get their times
good_spike_times = spike_times(spike_indices);
% and their clusters
good_spike_clusters = spike_clusters(spike_indices);
% and their probes!
good_probes = probes(good_cluster_indices);

% get only the first probe data 
PROBE_1_MAX = 551;
probe1_clusters_indices = find(good_spike_clusters > PROBE_1_MAX);
probe1_clusters = good_spike_clusters(probe1_clusters_indices);
probe1_spike_times = good_spike_times(probe1_clusters_indices);
%{
Structure:
    TODO: raster plot of spikes 
    clusters.       A collection of spikes to be analyzed together. 
                    They are considered to arise from a single neuron except where annotated otherwise.
       
        _phy_annotation.npy:  >= 2 is neuron to consider. 
            
    spikes.         Each element corresponds to one detected action potential
        
    
%}
%binaryBinnedSpikeCell = cellfun(@(x) arrayfun(@(a) a > 0, x, 'UniformOutput',false), binnedSpikeTimesCell, 'UniformOutput',false)


