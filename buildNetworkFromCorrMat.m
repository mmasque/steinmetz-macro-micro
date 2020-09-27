function network = buildNetworkFromCorrMat(corrMat, stdevTresh)
    % stdevTresh is how many stdeviations away from the mean to set the threshold to. 
    stdev = std(corrMat(:));
    thresh = stdevTresh * stdev;
    
    network = corrMat > thresh;
end