function indices = getClustersByDepthRange(low, hi, depths, indices)
    %{
        selects a subset of cluster indices by a range of depth
        given by low-hi (inclusive). 
        
        returns the ids of the selected clusters sorted by depth, with
        deepest (lowest depth value (quirk of the dataset)) first. 
    %}
    logicalRange = depths >= low & depths <= hi;
    depthsRange = depths(logicalRange);
    
    
    [~,sortIdx] = sort(depthsRange);
    % sort indices using the sorting index of depths
    indices = indices(sortIdx);
    
    
end