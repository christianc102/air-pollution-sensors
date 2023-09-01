function [RMSE,MPE] = rmse_mpe(recon,actual,mask)
% Calculates Root Mean Squared Error and 
% Mean Percentage Error of reconstruction with actual data
% recon is 213804 x 1, actual is 871 x 413 with NaN values

    % Apply mask on actual and reshape into vector
    X = zeros(length(mask(mask==1)),1);
    X = actual(mask==1);
    
    % No modifications needed for recon
    Y = recon; 
    
    % Compute RMSE 
    RMSE = sqrt(mean((Y(:)-X(:)).^2));
    
    % Compute MPE 
    accuracy = 100 * abs(Y-X) ./ X;
    MPE = mean(accuracy(accuracy < Inf));
end

