function inverseSquaredMat = inverseSquareProportion(corrMat, indices, distances)
    %not tested!
    inverseSquaredMat = corrMat(indices, indices) ./ (getMatOfDistances(distances, indices, 0)) .^ 2;
end

