function [sweeps,scrubbedtimes,times,scrubbedidx,nspk_og] = loadspikes_hpc(bindatdir,thischan,nSamples,tINFO,minrng,swprng)
% This function replaces a block of code in the original
% binarySpikeProcessing_HPC, in order to paralellize/streamline the code.
% The function takes as INPUTS:
% __bindatdir: path to directory containing binary spike files
% __channel: channel number to load
% __nSamples: number of samples in each waveform (e.g. 97)
% __tINFO: experiment start time (unix time format)
% __minrng: sample point range in which to look for minimum
% __swprng: sample point range to clip around spike minimum after alignment
%
% The function proceeds in the following way.
% 1. Load and reshape waveform array
% 2. Load spike timestamps array
% 3. Align spike timestamps to 7:30 am on first recording day
% 4. Eliminate noise (large waveform traces)
% 5. Align spikes to minimum
%
% OUTPUTS:
% __sweeps: waveform array without noise and with waveforms aligned to
%           minimum and clipped
% __scrubbedtimes: indices of timestamps disregarding noise traces
% __times: indices of all timestamps (including noise)
% __scrubbedidx: indices of waveforms in cleaned array
% __nspk_og: original number of waveforms recorded

fprintf('Loading spikes for channel %u.\n',thischan)
spkfname = dir([bindatdir filesep '*_Channel_' num2str(thischan) '_spikes*']);
[fid, ~] = fopen( [bindatdir filesep spkfname.name],'r');
spkfsize = spkfname.bytes;
if ~isequal(spkfsize/4,round(spkfsize/4))
    warning('WARNING: spike file size inconsistent with single-precision floats.');
end

% if the spikes file is very large, read in two chunks
if spkfsize > 2.5e10
    spikes = zeros(nSamples,spkfsize/(nSamples*4));
    fprintf('file too big! Loading in 2 chunks.\n');
    half_bytes = spkfsize/2; % half number of bytes in file
    for half = 1:2
        fprintf('chunk %u\n',half);
        offset_bytes = floor((half-1)*half_bytes);
        A_start = floor((half-1)*size(spikes,2)/2 + 1);
        A_stop  = floor(half*size(spikes,2)/2);
        fseek(fid, offset_bytes, 'bof');
        tempA = fread(fid, [nSamples (A_stop-A_start+1)], 'single');
        
        spikes( : , A_start : A_stop ) = tempA;
        clear tempA
    end
    disp('Done!');
else
    spikes = fread(fid, 'single=>single');
    spikes = reshape(spikes,nSamples,[]);
    fclose(fid);
    disp('Done!');
end

if ~isempty(spikes)
    
    % load spike timestamps
    fprintf('Loading spike timestamps for channel %u.\n',thischan)
    tmfname = dir([bindatdir filesep '*_Channel_' num2str(thischan) '_times*']);
    [fid, ~] = fopen( [bindatdir filesep tmfname.name],'r');
    B = [];
    B = fread(fid,'double');
    fclose(fid);
    times = B;
    % align to experiment start time
    toffset = mod(tINFO,24*3600) - (7.5*3600);
    times = times - tINFO + toffset; % adjust all spike times to account for offset from 7:30.
    clear B
    disp('Done!');
    
    % index and skip massive events
    disp('Scrubbing the data to remove obvious outliers. ');
    sz              = size( spikes, 1 );
    origindx        = [1:length(times)]';
    [~,poscol]      = find(spikes>400 | spikes <-600); % Filter out huge traces;
    poscol          = unique(poscol);
    scrubbedidx     = setdiff(origindx,poscol);
    nspk_og         = size(origindx); % HERE
    scrubbedspikes  = spikes(:,scrubbedidx);
    clear spikes poscol origindx
    scrubbedtimes   = times(scrubbedidx);
    [nrow,ncol] = size(scrubbedspikes);
    if ncol ~= sz && nrow == sz
        scrubbedspikes = scrubbedspikes';
    end
    disp('Outliers removed.');
    
    
    % Do a first alignment to minimum before interpolation.
    
    % center spikes on the minimum
    [~,minidx]  = min(scrubbedspikes(:,minrng),[],2);
    minidx      = minidx + (minrng(1) - 1);
    % align
    sweeps = alignrows(scrubbedspikes, minidx, swprng, nSamples - 1);
    clear scrubbedspikes
else
    disp('No spikes for this channel :-(');
    sweeps          = [];
    scrubbedtimes   = [];
    times           = [];
	scrubbedidx 	= [];
    nspk_og         = [];
end
