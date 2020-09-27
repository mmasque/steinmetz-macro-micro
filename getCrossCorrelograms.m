function crossCorrelograms = getCrossCorrelograms(binarisedMatrix, tauRange)
    crossCorrelograms = zeros(size(binarisedMatrix, 2),size(binarisedMatrix, 2), 2*tauRange+1);
    for i = 1:size(binarisedMatrix, 2)
        for j = 1:size(binarisedMatrix, 2)
            crossCorrelograms(i,j, :) = getCrossCorrelogram(binarisedMatrix(:,i), binarisedMatrix(:,j), tauRange);
        end
    end
    
end