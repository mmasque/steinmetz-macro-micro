%% constants
PATH = "allData/Cori_2016-12-14";
BINWIDTH = 0.005;
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
binnedSpikeTimesCell_5ms = spikeTimeBinning(spikeTimesCell, BINWIDTH);
%% get binary bins (spike or no spike within window)
binarySpikeTimesCell_5ms = spikeTimesToBinary(binnedSpikeTimesCell_5ms);
%% in matrix form
binarySpikeTimesMatrix = cell2mat(binnedSpikeTimesCell')';
%% get the correlation matrix
corrCoeffMatrix = getBinSpikesCorrCoeffs(binnedSpikeTimesCell_5ms(probeIndices));
%%
noNanCorrCoeffMatrix = corrCoeffMatrix(probeIndices, probeIndices);
%% get neuron group to analyse by depth
depthRangeProbeIndices = getClustersByDepthRange(0,3.7535e+03,probeDepths, probeIndices);
% depthRangeCorrCoeffMatrix = corrCoeffMatrix(depthRangeProbeIndices, depthRangeProbeIndices);
%% get frequency of simultaneous spikes
simSpikesMat = allSimSpikes(binnedSpikeTimesCell_5ms(probeIndices));
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
conns = goodConnections(1:9, :);
nValues = 2;
nStates = 2;
N_SAMPLES = 1000;
TAU = 1;
PATH = "tpms/TPM_CORI_2016_12_14_PROBE1_DEPSHUFFLE_TESTgoodConnections_randsample_10msbin_10x_diffsamples";
REPETITIONS = 10;
mkdir(PATH);
%tpms = zeros(size(conns,2)^nValues, size(conns,2)^nValues, size(conns,1), REPETITIONS);
%state_obs = zeros(size(conns,2)^nValues, size(conns,1), REPETITIONS);
for samp = 200:200:1000
for j = 1:REPETITIONS
        for i = 1:length(conns)
            tpm_data = [shuffleInTime(binarySpikeTimesCellProbeIndices_10ms{goodConnections(i,1)}, true); shuffleInTime(binarySpikeTimesCellProbeIndices_10ms{(goodConnections(i,2))}, true)]';
            [tpm, state_counters] = build_tpm_equalStates_randomsampling(tpm_data, TAU, nValues, samp);
            %tpms(:,:,i, j) = tpm;
            neurons = [goodConnections(i,1), goodConnections(i,2)];
            saveto = fullfile(PATH,strcat("TPM_n", string(goodConnections(i,1)), "_n", string(goodConnections(i,2)), "_samp", string(samp), "_rep", string(j)));

            save(saveto,"tpm", "state_counters", "nValues", "neurons");
            
            %state_obs(:, i, j) = countOccurrences(tpm_data, nValues, nStates);
        end
end
end
%% after running IIT in python, recover the phi values from results/split
PATH = "results/split/TPM_CORI_2016_12_14_PROBE1_goodConnections_randsample_10msbin_10x_diffsamples";
conns = goodConnections(1:9, :);
good_phis_10 = zeros(length(conns), REPETITIONS, 5);
ticklabels_good = strings(1,length(conns));
for k = 1:5
for j = 1:10
for i = 1:length(conns)
    fname = fullfile(PATH,strcat("TPM_n", string(conns(i,1)), "_n", string(conns(i,2)), "_samp",string(k*200), "_rep", string(j), ".mat"));
    load(fname);
    good_phis_10(i,j, k) = phi.phi;
    ticklabels_good(i) = strcat(string(conns(i,1)), ", ", string(conns(i,2)));
end
end
end
PATH = "results/split/TPM_CORI_2016_12_14_PROBE1_DEPSHUFFLE_TESTgoodConnections_randsample_10msbin_10x_diffsamples";
conns = goodConnections(1:9, :);
random_phis_10 = zeros(length(conns), REPETITIONS, 5);
ticklabels_random = strings(1,length(conns));
for k = 1:5
for j = 1:10
for i = 1:length(conns)
    fname = fullfile(PATH,strcat("TPM_n", string(conns(i,1)), "_n", string(conns(i,2)),"_samp",string(k*200), "_rep", string(j), ".mat"));
    load(fname);
    random_phis_10(i, j, k) = phi.phi;
    ticklabels_random(i) = strcat(string(conns(i,1)), ", ", string(conns(i,2)));
end
end
end
%%

figure;
s1 = subplot(1,2,1); 
   
    hold on
    errorbar(squeeze(mean(good_phis_5(:, :, 5),2)), squeeze(std(good_phis_5(:, :, 5),0,2)), 'o', 'LineWidth', 2);
    errorbar(squeeze(mean(good_phis_10(:, :, 5),2)), squeeze(std(good_phis_10(:, :, 5),0,2)), 'o', 'LineWidth', 2);
    hold off
    %hold on
    %boxplot(good_phis_10(:, :, 5)')
    %errorbar(squeeze(mean(good_phis_5(:, :, 5),2)), squeeze(std(good_phis_5(:, :, 5),0,2)), 'o', 'LineWidth', 2);

    xticks(1:length(good_phis));
    xticklabels(ticklabels_good);
    xtickangle(90);
    xlabel("Neuron pairs")
    ylabel("{\Phi}")
    set(gca,'FontSize',15)
    title({'Connected Neurons', 'randomly selected 1000 transitions/state for TPM', '10 iterations'});

s2 = subplot(1,2,2); 

    hold on
    errorbar(squeeze(mean(random_phis_5(:, :, 5),2)), squeeze(std(random_phis_5(:, :, 5),0,2)), 'o', 'LineWidth', 2);
    errorbar(squeeze(mean(random_phis_10(:, :, 5),2)), squeeze(std(random_phis_10(:, :, 5),0,2)), 'o', 'LineWidth', 2);
    hold off
    %boxplot(random_phis_10(:, :, 5)')
    %errorbar(squeeze(mean(random_phis_5(:, :, 5),2)), squeeze(std(random_phis_5(:, :, 5),0,2)), 'o', 'LineWidth', 2);

    xticks(1:length(random_phis));
    xticklabels(ticklabels_random);
    xtickangle(90);
    xlabel("Neuron pairs")
    ylabel("{\Phi}")
    set(gca,'FontSize',15)
    title({'Independently Time shuffled Conn Neurons', 'randomly selected 1000 transitions for TPM', '10 iterations'});

%linkaxes([s1, s2], 'y');
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
    if pairDistance < maxDistance && ~ismember([connectedNeurons(i, 1) connectedNeurons(i, 2)], goodConnections, "rows") && ~(connectedNeurons(i, 1) == connectedNeurons(i, 2))
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