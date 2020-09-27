function [network, indices] = findCrossCorrelogramsWithPeak(crossCorrelograms, thresh)
    network = zeros(size(crossCorrelograms,1), size(crossCorrelograms,1));
    indices = [];
    for i = 1:size(crossCorrelograms,1)
        for j = 1:size(crossCorrelograms,2)
            network(i,j) = sum(isoutlier(crossCorrelograms(i,j,:),'mean','ThresholdFactor',thresh)) > 0;
            if network(i,j)
                indices = [indices; [i,j]];
            end
        end
    end
end