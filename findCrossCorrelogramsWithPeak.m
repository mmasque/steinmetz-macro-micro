function [indices] = findCrossCorrelogramsWithPeak(crossCorrelograms, thresh)
    network = zeros(size(crossCorrelograms,1), size(crossCorrelograms,1));
    indices = [];
    for j = 1:size(crossCorrelograms,2)
        outliers = sum(isoutlier(crossCorrelograms(:, j),'mean','ThresholdFactor',thresh)) > 0;
        if outliers
            indices = [indices; j];
        end
    end
end