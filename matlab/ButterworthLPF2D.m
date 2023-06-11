function [Y,Yinv] = ButterworthLPF2D(X, D0, n)
%UNTITLED2 Summary of this function goes here
%   This cute little function generates a butterowrth filter based upon
%   input parameters D0 is the diameter, n is the order of the filter.
%   The size of the mask determined by the input image

    [M, N] = size(X);
    Y = zeros(M, N);
    [fy, fx] = ndgrid(0:M/2, 0:N/2);
    
    % Butterworth mask same center as Gaussian
    % n = 1 % filter order
    % Calculating Euclidean Distance
    D = sqrt(fy.^2 + fx.^2);
    Y = 1./(1 + (D./D0).^(2*n));
    % determining the filtering mask
    % Do symmetries
    Y(1:M/2+1, N:-1:N/2+2) = Y(1:M/2+1, 2:N/2);
    Y(M:-1:M/2+2, :) = Y(2:M/2, :);
    Yinv = 1 - Y;

end

