function distancesMat = getMatOfDistances(distancesVec, indicesVec, proximityLimit)
   % not tested!
    distancesMat = NaN(length(indicesVec), length(indicesVec));
    for i = 1:length(indicesVec)
        for j = 1:length(indicesVec)
            distance = abs(distancesVec(indicesVec(i)) - distancesVec(indicesVec(j)));
            if distance < proximityLimit
                distancesMat(i,j) = Inf;
            else
                distancesMat(i,j) = distance;
            end
        end
    end
end