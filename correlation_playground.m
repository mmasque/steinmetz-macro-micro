%% constants
PATH = "allData/Cori_2016-12-14";
BINWIDTH = 0.0005;
PROBE = 1;

%% data
probes = readNPY(strcat(PATH, "/clusters.probes.npy"));
depths = readNPY(strcat(PATH, "/clusters.depths.npy"));
%% get good indices per _phy.annotation label. 
[indices, clusters, times] = getGoodIndices(PATH);%% get data from needed probe
probeIndices = getIndicesForProbe(PROBE, probes, indices);
probeClusterIndices = find(ismember(clusters, probeIndices));
probeClusters = clusters(probeClusterIndices);
probeTimes = times(probeClusterIndices);
probeDepths = depths(probeIndices);
%% reformat the spikes and clusters data
spikeTimesCell = reformatSpikeTimes(probeIndices, probeTimes, probeClusters);
binnedSpikeTimesCell = spikeTimeBinning(spikeTimesCell, BINWIDTH);
%% get binary bins (spike or no spike within window)
binarySpikeTimesCell = spikeTimesToBinary(binnedSpikeTimesCell);
%% in matrix form
binarySpikeTimesMatrix = cell2mat(binnedSpikeTimesCell')';
%% get the correlation matrix
corrCoeffMatrix = getBinSpikesCorrCoeffs(binnedSpikeTimesCell);
%%
noNanCorrCoeffMatrix = corrCoeffMatrix(probeIndices, probeIndices);
%% get neuron group to analyse by depth
depthRangeProbeIndices = getClustersByDepthRange(0,3.7535e+03,probeDepths, probeIndices);
% depthRangeCorrCoeffMatrix = corrCoeffMatrix(depthRangeProbeIndices, depthRangeProbeIndices);

%% get cross correlograms
crossCorrelogram_1004_1027_0_00005 = getCrossCorrelogram(binnedSpikeTimesCell(:, 1027), binnedSpikeTimesCell(:, 1027), 50);
%%
%for i = 10:13
%for j = 10:13
for ind = 1:length(ccindices)
    i = ccindices(ind, 1);
    j = ccindices(ind,2);
    figure;
    bar(squeeze(crossCorrelograms(i,j,:)))
    title({"Cross Correlogram",strcat("neuron ", num2str(i), "(t) ", "and neuron ", num2str(j), "(t+d)"),strcat("Binwidth(ms): ", num2str(BINWIDTH * 1000))})
    xticks(1:5:41);
    xticklabels((-20:5:20) * BINWIDTH * 1000);
    xlabel("d (ms)")
    ylabel("counts")
    saveas(gcf, strcat("crossCorrelogramSelectedPairs_", num2str(i), "_", num2str(j), "_binwidth_ms_", num2str(BINWIDTH*1000), '.png'));
    close(gcf);
%end
end

