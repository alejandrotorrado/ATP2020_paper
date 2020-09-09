function [LD_data] = get_LD_times(tank,anim_ID)

fprintf('\nExtracting L/D transition times for animal: %s\n\n',anim_ID);
t00 = tic;

% find the relevant directories in the data tank
tank_dirs = dir([tank filesep 'Chronic-*']);

% initialize variables
ON_TIMES = [];
OFF_TIMES = [];
block_info = [];
block_dirs = {};

% loop through tank directories
for ii = 1:numel(tank_dirs)
    
    % find the corresponding block directory based on tank name
    block_dir = [tank filesep tank_dirs(ii).name];
    block_dirs{ii} = block_dir;
    
    fprintf('\n  Block %u, folder: %s\n',ii,block_dir);
    
    % Use the TDT function TDTbin2mat to extract the data from the TEV
    % files
    t0 = tic;
    data_this_block = TDTbin2mat(block_dir,'TYPE',{'epocs'},'STORE','Lite','VERBOSE',0);
    t1 = toc(t0);
    fprintf('  Block %u took %.2f seconds.\n',ii,t1);
    
    % if there are L-D transition times stored here
    if isfield(data_this_block.epocs,'Lite')
        
        % attempt to get them
        try
            tmp_lite = data_this_block.epocs.Lite;
        catch
            keyboard;
        end
        
        % onset are the actual times, data is the voltage trace
        tmp_onoff = tmp_lite.data;
        tmp_times = tmp_lite.onset;
        
        % get TSQ timestamp for start of this recording session
        t_start_tmp = checkTSQ(block_dir);
        block_start_time = unixtime(datevec(TimezoneConvert(datestr(unixtime(t_start_tmp))...
            ,'UTC','America/New_York')));
        
        % adjust L-D times based on that
        tmp_ON_TIMES = tmp_times(tmp_onoff==1) + block_start_time;
        tmp_OFF_TIMES = tmp_times(tmp_onoff==0) + block_start_time;
        
        % store
        ON_TIMES = [ON_TIMES; tmp_ON_TIMES];
        OFF_TIMES = [OFF_TIMES; tmp_OFF_TIMES];
        block_info = [block_info; [ii block_start_time]];
       
    else
         fprintf('No Lite data for this block.\n\n');
    end
    
end

% create output variables
LD_data.ON_TIMES    = ON_TIMES;
LD_data.OFF_TIMES   = OFF_TIMES;    
LD_data.animal      = anim_ID;
LD_data.block_info  = block_info;
LD_data.block_dirs  = block_dirs;

t01 = toc(t00);
fprintf('Total time elapsed: %.2f seconds.\n',t01);



