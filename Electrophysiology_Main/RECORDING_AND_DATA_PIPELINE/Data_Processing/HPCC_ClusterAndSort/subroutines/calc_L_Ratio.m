function [outLR, mdist] = calc_L_Ratio(features, clustSpikes)
% output = L_Ratio(FD, ClusterSpikes)
%
% L-ratio
% Measure of cluster quality
%
% Inputs:   FD:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into FD which lists spikes from the cell whose quality is to be evaluated.
%
% Output: a structure containing three components
%           Lratio, L, df

% find # of spikes in this cluster
[nSpikes, nFeat] = size(features);

nClusterSpikes = length(clustSpikes);

% mark spikes which are not cluster members
NoiseSpikes = setdiff(1:nSpikes, clustSpikes);

mdist = mahal(features, features(clustSpikes,:));
df = size(features,2);

L = sum(1-chi2cdf(mdist(NoiseSpikes),df));
Lratio = L/nClusterSpikes;

outLR.L = L;
outLR.Lratio = Lratio;
outLR.df = nFeat;