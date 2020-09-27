function network = getNetworkFromCSV(csvFileName, dim)
    network = zeros(dim,dim);
    connectionData = readmatrix(csvFileName);
    for i = 1:length(connectionData)
        if connectionData(i,4) == 1
            network(connectionData(i,1), connectionData(i,2)) = 1;
        end
    end
    
end