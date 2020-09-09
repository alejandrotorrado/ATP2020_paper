function [t0] = checkTSQ(direc)
% 5_17_16 This script will check each block folder for the presence of a
% TSQ file. The TSQ file, among other things, contains the unix time stamp
% of the start of the recording. Therefor, this file will only exist in the
% first block folder of a recording, or the first blcok folder following a
% crash. To figure out subsequent block start times, we've added a time
% counter in channels 65-68. This counter starts at 0 at the beginning of
% the recording, and thus, in combination with the TSQ file, can be used to
% make accurate time stamps for the beginning of each block.
%
% This script will check each folder for a TSQ file, and, if it is present,
% return the start time (unix time) of the recording. If there is no TSQ
% file, this will return a 0.

tsq = dir(fullfile(direc, '*.tsq'));

if ~isempty(tsq);
    
    
    EVMARK_STARTBLOCK   = hex2dec('0001');
    EVMARK_STOPBLOCK    = hex2dec('0002');
    
    tsqdat = fopen([direc filesep tsq.name], 'rb');
    fseek(tsqdat, 48, 'bof');  code1 = fread(tsqdat, 1, '*int32');
    assert(code1 == EVMARK_STARTBLOCK);
    fseek(tsqdat, 56, 'bof');
    start_time = fread(tsqdat, 1, '*double');
    
    fclose(tsqdat);
    
    %
    %     % read stop time
    %     fseek(tsqdat, -32, 'eof'); code2 = fread(tsqdat, 1, '*int32');
    %     assert(code2 == EVMARK_STOPBLOCK);
    %     fseek(tsqdat, -24, 'eof'); stop_time = fread(tsqdat, 1, '*double');
    %     % total duration for data size estimation
    %     duration = stop_time-start_time;
    %
    t0 = start_time;
else
    t0 = 0;
end