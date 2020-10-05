%% constants
PATH = "allData/Cori_2016-12-14";
BINWIDTH = 0.001;
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
binnedSpikeTimesCell_10ms = spikeTimeBinning(spikeTimesCell, 0.01);
%% get binary bins (spike or no spike within window)
binarySpikeTimesCell_10ms = spikeTimesToBinary(binnedSpikeTimesCell_10ms);
%% in matrix form
binarySpikeTimesMatrix = cell2mat(binnedSpikeTimesCell')';
%% get the correlation matrix
corrCoeffMatrix = getBinSpikesCorrCoeffs(binnedSpikeTimesCell_10ms(probeIndices));
%%
noNanCorrCoeffMatrix = corrCoeffMatrix(probeIndices, probeIndices);
%% get neuron group to analyse by depth
depthRangeProbeIndices = getClustersByDepthRange(0,3.7535e+03,probeDepths, probeIndices);
% depthRangeCorrCoeffMatrix = corrCoeffMatrix(depthRangeProbeIndices, depthRangeProbeIndices);
%% get frequency of simultaneous spikes
simSpikesMat = allSimSpikes(binnedSpikeTimesCell_10ms(probeIndices));
%% get the binarySpikeTimesCell with only probeIndices
binnedSpikeTimesCellProbeIndices = binnedSpikeTimesCell(probeIndices);
%% binarise the probe indices
binarySpikeTimesCellProbeIndices = spikeTimesToBinary(binnedSpikeTimesCellProbeIndices);
%% get indices of potentially coupled neurons 
network = buildNetworkFromCorrMat(simSpikesMat, 4);
imagesc(network)
connectedNeurons = getConnectedNeurons(network);
%% get cross correlograms of the potentially coupled neurons
% these cross correlograms are indexed by the pairs of connectedNeurons
TAURANGE = 15;
crossCorrelograms = zeros(2*TAURANGE + 1,length(connectedNeurons));
for i = 1:length(connectedNeurons)
    crossCorrelograms(:, i) = getCrossCorrelogram(binarySpikeTimesCellProbeIndices{connectedNeurons(i, 1)},binarySpikeTimesCellProbeIndices{connectedNeurons(i, 2)}, TAURANGE);
end
%% filter the cross correlograms by finding if there is a peak
connectedNeuronsFilteredByCC = connectedNeurons(findCrossCorrelogramsWithPeak(crossCorrelograms, 3.5), :);
% remove self connections 
connectedNeuronsFilteredByCC = connectedNeuronsFilteredByCC(connectedNeuronsFilteredByCC(:, 1) ~= connectedNeuronsFilteredByCC(:, 2), :);
ccindices = find(ismember(connectedNeurons,connectedNeuronsFilteredByCC, 'rows'));
%% plotting cross-correlograms
%for i = 10:13
%for j = 10:13
PATH = 'figures/crossCorrelograms_1msbin/';
mkdir(PATH);
for ind = ccindices'
    
    figure;
    bar(squeeze(crossCorrelograms(:, ind)))
    title({"Cross Correlogram",strcat("neuron ", num2str(connectedNeurons(ind, 1)), "(t) ", "and neuron ", num2str(connectedNeurons(ind, 2)), "(t+d)"),strcat("Binwidth(ms): ", num2str(BINWIDTH * 1000))})
    xticks(1:5:31);
    xticklabels((-15:5:15) * (BINWIDTH * 1000));
    xlabel("d (ms)")
    ylabel("counts")
    saveas(gcf, strcat(PATH,"crossCorrelogramSelectedPairs_", num2str(connectedNeurons(ind, 1)), "_", num2str(connectedNeurons(ind, 2)), "_binwidth_ms_", num2str(BINWIDTH*1000), '.png'));
    close(gcf);
%end
end

%% get TPMs of connected neurons
conns = randomNeuronPairs(1:9, :);
nValues = 2;
nStates = 2;
N_SAMPLES = 500;
TAU = 1;
PATH = "tpms/TPM_CORI_2016_12_14_PROBE1_randomNeuronPairs_randsample_1msbin";
mkdir(PATH);
tpms = zeros(size(conns,2)^nValues, size(conns,2)^nValues, size(conns,1));
state_obs = zeros(size(conns,2)^nValues, size(conns,1));
for i = 1:length(conns)
    tpm_data = [binarySpikeTimesCellProbeIndices{conns(i,1)}; binarySpikeTimesCellProbeIndices{conns(i,2)}]';
    [tpm, state_counters] = build_tpm_equalStates_randomsampling(tpm_data, TAU, nValues, N_SAMPLES);
    tpms(:,:,i) = tpm;
    neurons = [conns(i,1), conns(i,2)];
    saveto = fullfile(PATH,strcat("TPM_n", string(conns(i,1)), "_n", string(conns(i,2))));
    
    save(saveto,"tpm", "state_counters", "nValues", "neurons");
    
    state_obs(:, i) = countOccurrences(tpm_data, nValues, nStates);
end

%% after running IIT in python, recover the phi values from results/split
PATH = "results/split/TPM_CORI_2016_12_14_PROBE1_randomNeuronPairs_randsample_1msbin";
conns = randomNeuronPairs(1:9, :);
phis = zeros(1,length(conns)*2);
ticklabels = strings(1,length(conns)*2);
for i = 1:length(conns)
    fname = fullfile(PATH,strcat("TPM_n", string(conns(i,1)), "_n", string(conns(i,2))));
    load(fname);
    phis(i) = phi.phi;
    ticklabels(i) = strcat(string(conns(i,1)), ", ", string(conns(i,2)));
end

PATH = "results/split/TPM_CORI_2016_12_14_PROBE1_Manual_network_pairs_randsample_1msbin";
conns = goodConnections;
for i = 1:length(conns)
    fname = fullfile(PATH,strcat("TPM_n", string(conns(i,1)), "_n", string(conns(i,2))));
    load(fname);
    phis(i + length(conns)) = phi.phi;
    ticklabels(i + length(conns)) = strcat(string(conns(i,1)), ", ", string(conns(i,2)));
end

colors = zeros(length(phis), 3);
colors(1:length(conns), :) = ones(length(conns), 3) .*  [1 0 0];
colors(1+length(conns):length(phis), :) = ones(length(conns), 3) .*[0 0 1];
figure;
scatter(1:length(phis), phis, 80, colors, 'filled');
xticks(1:length(phis));
xticklabels(ticklabels);
xtickangle(90);
xlabel("Neuron pairs")
ylabel("{\Phi}")
set(gca,'FontSize',15)
title("Phi values of connected neuron pairs vs. randomly selected neuron pairs");
% add correct legend label (hacky)
hold on; scatter([], [], 'blue', 'filled');hold off
legend(["Random pairs", "Connected neurons"]);

%% get some neuron pairs that are at a similar distance to the mean distance between connected pairs. 

% compute mean distance between connected pairs
totalDistance = 0;
maxDistance = 0;
for i = 1:length(goodConnections)
    distance = abs(probeDepths(goodConnections(i, 1)) - probeDepths(goodConnections(i, 2)));
    maxDistance = max(maxDistance, distance);
    totalDistance = totalDistance + distance;
end
meanDistance = totalDistance/length(goodConnections);

% find neuron pairs with high coupled firing rates that have similar
% distances to the mean. 
randomNeuronPairs = [];
rng(0,'twister');   % repeatable randoms
r = randi([1 length(connectedNeurons)],1,length(connectedNeurons));
rindex = 1;
while size(randomNeuronPairs, 1) < length(goodConnections)
    
    i = r(rindex);
    pairDistance = abs(probeDepths(connectedNeurons(i, 1)) - probeDepths(connectedNeurons(i, 2)));
    if pairDistance < meanDistance + maxDistance && ~ismember([connectedNeurons(i, 1) connectedNeurons(i, 2)], goodConnections, "rows") && ~(connectedNeurons(i, 1) == connectedNeurons(i, 2))
        randomNeuronPairs = [randomNeuronPairs; [connectedNeurons(i, 1) connectedNeurons(i, 2)]];
    end
    rindex = rindex + 1;
end

%% compute average TPM
path = "tpms/TPM_CORI_2016_12_14_PROBE1_Manual_network_pairs_randsample";
sumTPM = zeros(4,4);
for i = 1:length(goodConnections)
    fname = strcat(path, "/TPM_n", string(goodConnections(i,1)), "_n", string(goodConnections(i,2)), ".mat");
    tpm = load(fname).tpm;
    sumTPM = sumTPM + tpm;
end
meanTPMGoodConnections = sumTPM./length(goodConnections);

%% plot a tpm, from the variable tpm, for neuron pair 'neurons', binwidth BINWIDTH
TITLE = {"Cori 2016 12 14 Probe 1","Transition probability matrix of pair of neurons, mean across identified connected pairs",strcat("binwidth ", string(BINWIDTH*1000), "ms, random sampling")};
h = heatmap(meanTPMGoodConnections);
h.Title = TITLE;
h.XData = XData;
h.YData = YData;
h.ColorLimits = [0 1];
h.GridVisible = 0;
h.YLabel = "states of Neuron A Neuron B at t-1";
h.XLabel = "states of Neuron A Neuron B at t";
h.FontSize = 14;


%% compute phi for subgraphs
SAMPLES = 500;
nVals = 2;
TAU = 1;
subgraphs = {[10 323], [38 307], [291 232 171], [348, 55, 305], [168,143,331,293]};
tpms_10ms = cell(length(subgraphs),1);
PATH = "tpms/TPM_CORI_2016_12_14_PROBE1_Manual_network_groups_10ms";
mkdir(PATH);
for i =  1:length(subgraphs)
    neurons = subgraphs{i};
    nStates = length(sub);
    tpm_data = [];
    fname = "TPM";
    for node = neurons
        tpm_data = [tpm_data;binarySpikeTimesCellProbeIndices_10ms{node}];
        fname = strcat(fname, "_n", string(node));
    end
    fname = strcat(fname, ".mat");
    tpm_data = tpm_data';
    [tpm, state_counters] = build_tpm_equalStates_randomsampling(tpm_data, TAU, nVals, SAMPLES);
    tpms_10ms{i} = tpm;
    save(fullfile(PATH, fname),"tpm", "state_counters", "nValues", "neurons");

end