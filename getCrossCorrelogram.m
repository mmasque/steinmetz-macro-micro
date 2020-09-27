function cross_correlogram = getCrossCorrelogram(a,b,tauRange)
    % from -tauRange to +tauRange
    cross_correlogram = zeros(1, 2*tauRange);
    for i = -1*tauRange:1:tauRange
        
        %try
            shiftedB = circshift(b, i);
            if i > 0
                shiftedB(1:i) = 0;%zeros(1:i);
            elseif i < 0
                shiftedB(length(shiftedB)+i+1:end) = 0;
            end
            cross_correlogram(i + tauRange + 1) = nnz(bitand(a, shiftedB)); %number of 1s in the bitwise and of both vectors.
        %catch
         %   ME = MException('CrossCorrelogram:tauRangeExceedsVector', ...
          %  "tauRange exceeds vector length.");
           % throw(ME);
        %end
    end
    % reverse it at the end, now at index t it is matching a(i) with b(i+t).
    cross_correlogram = flip(cross_correlogram);
end

