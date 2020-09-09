function IsoDist = Iso_Distance(features, clustSpikes)

% IsoDist = IsolationDistance(FD, ClusterSpikes)
%
% Isolation Distance
% Measure of cluster quality
%
% Inputs:   FD:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into FD which lists spikes from the cell whose quality is to be evaluated.
%
% Created by Ken Harris
% 
% Code by ADR 2012/12, from earlier versions
[nSpikes, nFeat] = size(features);

nClusterSpikes = length(clustSpikes);

if nClusterSpikes > nSpikes/2
    IsoDist = nan; % not enough out-of-cluster-spikes - IsoD undefined
    return
end

InClu = clustSpikes;
OutClu = setdiff(1:nSpikes, clustSpikes);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
m = mahal(features, features(clustSpikes,:));

mNoise = m(OutClu); % mahal dist of all other spikes

% calculate point where mD of other spikes = n of this cell
sorted = sort(mNoise);
IsoDist = sorted(nClusterSpikes);
