function [indices] = getIndicesForProbe(probeSelected, probesArray, goodIndices)

    goodProbes = probesArray(goodIndices);
    indices = goodIndices((goodProbes==probeSelected));
    
end