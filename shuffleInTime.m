function shuffled = shuffleInTime(dataset, together)
    % dataset -- channels x time
    % together -- shuffle channels together in time vs independently
    
    if together
        shuffled = dataset(:, randperm(size(dataset,2)));
    else
        shuffled = zeros(size(dataset));
        for i = 1:size(dataset,1)
            p = randperm(size(dataset,2));
            shuffled(i, :) = dataset(i, p);
        end
    end 
    
end
