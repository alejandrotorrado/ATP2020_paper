

function [Crik_data] = get_Crik_times(tank,anim_ID)

fprintf('\nExtracting Cam1 transition times for animal: %s\n\n',anim_ID);
t00 = tic;

% tank_dirs = dir([tank filesep 'Chronic-*']);
tank_dirs = dir([tank filesep 'Crickets_BJL-*']);

ON_TIMES = [];
OFF_TIMES = [];
block_info = [];
block_dirs = {};
blockcam_timestamps = {};

for ii = 1:numel(tank_dirs)
    
    block_dir = [tank filesep tank_dirs(ii).name];
    block_dirs{ii} = block_dir;
    
    fprintf('\n  Block %u, folder: %s\n',ii,block_dir);
    t0 = tic;
    data_this_block = TDTbin2mat(block_dir,'TYPE',1,'STORE','Crik','VERBOSE',1); % change Lite  
    %Previously used TYPE,{'epocs','snips'}
    t1 = toc(t0);
    fprintf('  Block %u took %.2f seconds.\n',ii,t1);
    
    if isfield(data_this_block.epocs,'Crik') % Change Lite
        
        
        try
            tmp_Crik = data_this_block.epocs.Crik; % change Lite
        catch
            keyboard;
        end
        
        
        
        tmp_onoff = tmp_Crik.data;
        tmp_times = tmp_Crik.onset;
        
        t_start_tmp = checkTSQ(block_dir);
        block_start_time = unixtime(datevec(TimezoneConvert(datetime(unixtime(t_start_tmp))...
            ,'UTC','America/New_York')));
        
        tmp_ON_TIMES = tmp_times(tmp_onoff==1) + block_start_time;
        tmp_OFF_TIMES = tmp_times(tmp_onoff==0) + block_start_time;
        
        ON_TIMES = [ON_TIMES; tmp_ON_TIMES];
        OFF_TIMES = [OFF_TIMES; tmp_OFF_TIMES];
        block_info = [block_info; [ii block_start_time]];
       
    else
         fprintf('No Cam1 data for this block.\n\n');
    end
    if exist('tmp_times')
    blockcam_timestamps{ii,1} = tmp_times
    end
end




Crik_data.ON_TIMES    = ON_TIMES;
Crik_data.OFF_TIMES   = OFF_TIMES;    
Crik_data.animal      = anim_ID;
Crik_data.block_info  = block_info;
Crik_data.block_dirs  = block_dirs;
% Crik_data.additional = data_this_block;
Crik_data.timestamps = blockcam_timestamps;
t01 = toc(t00);
fprintf('Total time elapsed: %.2f seconds.\n',t01);



