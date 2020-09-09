function [isi, varargout] = isidist (times,varargin)
% isidist.
% Interspike interval distribution calculation routine. Called by various 
% functions that process neural spiking data. Also used for plotting 
% features of neuronal activity. 
%
% INPUTS:
% times:  A list of all spike times occuring in a window of tims. This is 
%           the only required input. Spike times must be in seconds, but
%           can be either unix time or relative to an experiment start time
%           or time zone standard
%
%
% isi = isidist(times, PropName1, PropValue1.....)
%
% Valid field/value pairs are:
%       'edges'         Bin edges for histc binning of ISIs. Typically
%                       edges may be [0:0.001:1]
%
%       'cutoff'        The maximum time to be considered when returning
%                       the ISI distribution, binned data, and 
%                       contamination percentage.
%
% OUTPUTS:
% out1 = isidist(times) returns the interspike intervals of the list of
% times. Equivalent to "diff(times)". 
%
% [out1, out2] = isidist(...,'edges', [b0:n:b1]) also returns the histc
% binned data from out1.
%
% [out1, out2, out3] = isidist(times, ...) will calculate the percentage of
% total ISIs that are less than 3 milliseconds. See note below.
%
% KBH July 2016

edges   = [];
cutoff  = [];

for ee = 1:2:length(varargin);
    
    switch lower(varargin{ee})
        case 'edges'
            edges = varargin{ee+1};
        case 'cutoff'
           cutoff = varargin{ee+1};
    end
    
end

isi = diff(times);

% Only consider ISIs smaller than the cut-off value. 
if ~isempty(cutoff);
    isi = isi(isi<cutoff);
end

% Bin the ISI for easy plotting.
if ~isempty(edges);
    isiperbin = histc(isi,edges);  
else
    isiperbin = [];
end

% Calculate the probable error rate of the spikes in this cluster. Single
% units are unlikely to fire with a sub-three millisecond ISI, thus the
% percent of the distribution falling under that line can be used as an
% indication of error rate.
contamination = sum(isi<0.003)/numel(isi);


% - - - - - - - - - - - 

nout = nargout -1 ;

varlist = {'isiperbin','contamination'};

for ee = 1:nargout-1;
    
    varargout(ee) = {eval(varlist{ee})};
    
end



