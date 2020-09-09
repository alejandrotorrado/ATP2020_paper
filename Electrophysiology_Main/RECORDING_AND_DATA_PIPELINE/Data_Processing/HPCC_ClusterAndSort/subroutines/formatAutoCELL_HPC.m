function [newcell] = formatAutoCELL_HPC(celldat,nchunks,remove_badquals)
%
% formatAutoCELL (ATP - September 2016, as part of data pipeline
% consolidation effort).
%
% This function takes the output of the auto-sorting script and formats it
% in a way that is intuitive to use and compatible with our previous data
% structures.
%
% In essence, there are two possible cases:
%   1. the data was clustered over the whole experiment.
%      In this case, the data will be packaged as a CELL structure. Each
%      element in the structure (e.g. CELL(4)) will correspond to a cluster
%      identified by KlustaKwik, and will have the usual fields (quality,
%      channel, time, etc.). Thus the 'newcell' output variable will be a
%      list of all cells recorded (the usual CELL variable).
%
%   2. the data was clustered in chunks
%      In this case the output variable 'newcell' will take the form of a
%      structure that we will call CHUNKS. Each element (e.g. CHUNKS(1)) in
%      this will represent a chunk over which the clustering happened. It
%      will have fields indicating the start and end time of that chunk, as
%      well as a CELL field which will correspond to the usual CELL
%      variable (i.e. a structure-list of all identified clusters, with all
%      the relevant fields).

if nargin < 3
    remove_badquals = 1;
end


% create output structure
newcell = struct();
tempcell = [];

% if nchunks is 0, it means the dataset was clustered over the duration of
% the experiment, and the output structure should be a list of all cells
% identified during the clustering, named CELL.
if nchunks == 0
    % loop through channels
    for cc = 1:size(celldat,2);
        % if channel is not empty (i.e. not an EMG or ref channel)
        if ~isempty(celldat(cc).block)
            
            % store info about the channel the cell was recorded on
            thischan = celldat(cc).channel;
            
            % store info about all clusters identified
            allcells = celldat(cc).block.clust;
            
            for ii = 1:size(allcells,2)
                % all cells in this channel have same 'channel' field
                allcells(ii).channel = thischan;
                % store the p_save variable per cell
                allcells(ii).psort = celldat(cc).psave(ii,:);
            end
            
            % append list of clusters to the output structure
            if ~isfield(newcell,'quality')
                newcell = allcells;
            else
                newcell = [newcell allcells];
            end
        end
    end
    
    % remove quality 4 units
    if remove_badquals
        quals = [newcell.quality];
        kills = find(quals == 4);
        newcell(kills) = [];
    end
    
    % remove fields we don't use/need
    newcell = rmfield(newcell,{'dists','points','ctr'});
    
    
    
    % otherwise, there are several chunks of experiment. The ouptut structure
    % will be called CHUNKS and be structured such that CHUNKS(1) contains all
    % info about the first clustered piece. This will include the start and end
    % times of that piece, and all the info about all clusters identified by
    % KlustaKwik (i.e. a CELL structure).
elseif nchunks>0
    % loop through chunks
    for chunk = 1:nchunks
        
        % save chunk start and end time
        tmpblock = []; tmpcount = 0;
        while isempty(tmpblock)
            tmpcount = tmpcount + 1;
            try
            tmpblock = celldat(tmpcount).block{chunk};
            catch
                keyboard;
            end
        end
        newcell(chunk).START = celldat(tmpcount).block{chunk}.start;
        newcell(chunk).END = celldat(tmpcount).block{chunk}.end;
        
        % loop through channels
        for cc = 1:size(celldat,2);
            % if channel is not empty (i.e. not an EMG or ref channel)
            if ~isempty(celldat(cc).channel) && ~isempty(celldat(cc).block{chunk})
                
                % find channel
                thischan = celldat(cc).channel;
                
                % get all clusters for this chunk
				if ~isempty(celldat(cc).block{chunk})
                	try
                    	allcells = celldat(cc).block{chunk}.clust;
                	catch
						disp('Keyboard in formatAutoCell - chunks part');
                    	keyboard;
                	end
                
                	for ii = 1:size(allcells,2)
                    	% assign to each cluster the appropriate channel field
                    	allcells(ii).channel = thischan;
                    	% store the p_save variable per cell
                    	allcells(ii).psort = celldat(cc).psave{chunk}(ii,:);
                	end
                
                	% append to temporary structure-list variable
                	if ~isfield(tempcell,'quality')
                    	tempcell = allcells;
                	else
                   	 	tempcell = [tempcell allcells];
                	end
				else
					fprintf('Chunk %u is empty. Skipping.',chunk);
				end
            end
        end
        
        
        % remove quality 4 units
        if remove_badquals
            quals = [tempcell.quality];
            kills = find(quals == 4);
            tempcell(kills) = [];
        end
        
        % remove fields we don't use/need
        tempcell = rmfield(tempcell,{'dists','points','ctr'});
        
        % the CELL variable for each chunk is that temporary structure
        % containing the list of all clusters in that chunk
        newcell(chunk).CELL = tempcell;
        tempcell = [];
    end
end




