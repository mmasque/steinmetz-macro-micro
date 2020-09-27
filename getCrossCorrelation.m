function crossCorrelationMatrix = getCrossCorrelation(spikemat, maximumShift)
    crossCorrelationMatrix = NaN(size(spikemat, 2), size(spikemat, 2));
    subtractedMean = zeros(size(spikemat));
    for i = 1:size(spikemat, 2)
        subtractedMean(:, i) = spikemat(:,i) - mean(spikemat(:,i));
    end
    for i = 1:size(subtractedMean, 2)
        for ii = 1:size(subtractedMean,2)
            crossCorrelationMatrix(i,ii) = max(xcorr(subtractedMean(:, i), subtractedMean(:, ii), maximumShift, 'normalized'));
        end
    end
end