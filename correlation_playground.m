%% constants
PATH = "allData/Cori_2016-12-14";
BINWIDTH = 0.1;
PROBE = 1;

%% data
probes = readNPY(strcat(PATH, "/clusters.probes.npy"));
depths = readNPY(strcat(PATH, "/clusters.depths.npy"));
%% get good indices per _phy.annotation label. 
[indices, clusters, times] = getGoodIndices(PATH);
%% get data from needed probe
probeIndices = getIndicesForProbe(PROBE, probes, indices);
probeClusterIndices = find(ismember(clusters, probeIndices));
probeClusters = clusters(probeClusterIndices);
probeTimes = times(probeClusterIndices);
probeDepths = depths(probeIndices);
%% reformat the spikes and clusters data
spikeTimesCell = reformatSpikeTimes(probeIndices, probeTimes, probeClusters);
binnedSpikeTimesCell = spikeTimeBinning(spikeTimesCell, BINWIDTH);
%% get binary bins (spike or no spike within window)
binaryBinnedSpikeCell = cellfun(@(x) arrayfun(@(a) a > 0, x, 'UniformOutput',false), binnedSpikeTimesCell, 'UniformOutput',false);
%% get the correlation matrix
corrCoeffMatrix = getBinSpikesCorrCoeffs(binnedSpikeTimesCell);
%%
noNanCorrCoeffMatrix = corrCoeffMatrix(probeIndices, probeIndices);
%% get neuron group to analyse by depth
depthRangeProbeIndices = getClustersByDepthRange(0,300,probeDepths, probeIndices);
depthRangeCorrCoeffMatrix = corrCoeffMatrix(depthRangeProbeIndices, depthRangeProbeIndices);
