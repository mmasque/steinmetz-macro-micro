function simSpikes = countSimultaneousSpikes(a,b)
    simSpikes = 0;
    for i = 1:length(a)
        simSpikes = simSpikes + min(a(i), b(i));
    end
    
end
