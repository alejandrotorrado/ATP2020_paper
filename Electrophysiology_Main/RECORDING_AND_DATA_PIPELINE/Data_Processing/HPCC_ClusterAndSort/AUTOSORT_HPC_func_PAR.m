function AUTOSORT_HPC_func_PAR(animal,chanidx,nSamples,dephem,chanstart,binsize,startHere,endHere)

% autosort_hpc
%
% Clustering and auto-sorting script designed to run on Brandeis HPC
% cluster. This takes binary spikes/times files and returns a cleaned up
% CELL variable.
%
% ____________________________ IMPORTANT NOTE ____________________________ 
%
% The input arguments to this function are passed down from job files. An
% example job file can be found in:
%
% RECORDING_AND_DATA_PIPELINE/Data_Processing/HPCC_ClusterAndSort/ExampleJobs/
%
% ________________________________________________________________________
%
% This code uses KlustaKwik -> make sure this is available on the cluster.
% Save in a local directory and specify that location in
%
%
% FUNCTIONS AND SUBFCUNTIONS CALLED (ADD TO CLUSTER PATH):
%   - alignandcluster
%       - alignrows
%       - klustakwikRJ
%           - klustakwik_matlab
%   - alignrows
%   - autoClusterEval
%       - overclusteringcheck
%       - findpeaks
%       - isidist
%   - notchunking
%       - also needs alignrows, alignandcluster, autoClusterEval
%   - formatAutoCELL
%   - CHUNKS_addFields

% INPUTS - pass as arguments?
if binsize ~= 0
	disp('Autoclustering/sorting on chunked dataset.');
    startEdge = startHere;
    endEdge = endHere;
else
	disp('Autoclustering/sorting on continuous dataset.');
	startEdge = 0;
	endEdge = 0;
end



% sDir used to be /data/netapp/atorrpac on HPC64
% Code now updated for HPCC - 12/14/2018 ATP
sDir = ['/work/atorrpac/' animal filesep]; % folder containing BinarySpikes folder
addpath(sDir,'-end');
saveDir = ['/work/atorrpac/' animal filesep]; % output directory


% ADD PATH TO SUBROUTINES AND OTHER MATLAB FILES
% this location unchanged in HPCC - 12/14/2018 ATP
% This specifies the location of all the code needed to run the clustering
% and sorting. Change this to make it work for you. Note that you can have
% everything organized in sub-directories and just specify the top-most
% level folder here
CODE_PATH = '/home/atorrpac/Utilities/hpc_files'; 
addpath(genpath(CODE_PATH),'-end');


% call to the main script that does the calculations
binarySpikeProcessing_HPC_rForest_PAR(sDir,saveDir,animal,chanidx,nSamples,dephem,...
    chanstart,binsize,startEdge,endEdge);
