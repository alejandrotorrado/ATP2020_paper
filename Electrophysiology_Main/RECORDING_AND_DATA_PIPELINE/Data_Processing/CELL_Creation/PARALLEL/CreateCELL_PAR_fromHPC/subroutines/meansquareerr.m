function [mse] = meansquareerr (vec1, vec2)
% meansquareerr - is a simple calculation of the mean of the sum squared  
% residuals (deviations between corresponding points in vec1 and vec2).
% This is a commonly used estimate of the goodness of fit between a
% dataset (one of the inputs) and a model (the other input). 
%
% INPUTS:
%       vec1: a vector of N points contianing numeric data.
%       vec2: as numeric vector of the same length as in vec1
%
% OUTPUTS:
%       mse: the mean of the residual sum of squares. This is the error 
%           estimate.
%
% KBH July 2016
%

residuals = (vec1 - vec2).^2 ;

mse = mean(residuals);


