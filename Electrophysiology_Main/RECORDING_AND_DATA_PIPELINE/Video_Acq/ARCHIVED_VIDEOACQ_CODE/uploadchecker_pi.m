function uploadchecker_pi(timestamp,svdir,vrb)

% UPLOADCHECKER_v3: This version is for use with the RIGBOT. It will save a
% timestamp to the raspberry pi (path fed in via argument piDir) every time
% it is called by the main script. If the macmini file is found to have
% more than 1000 lines it will be wiped (except for the last timestamp).
% The timestamp argument should be a UNIX time.
% Set vrb = 1 for verbose mode.


if nargin < 3
    vrb = 0; % disable verbose mode as default
end

if ~exist(svdir,'dir')
    error('UPLOADCHECKER_PI: could not find RPi.');
else
    if vrb == 1, disp('UPLOADCHECKER_PI: found the RPi!'); end
end

if exist([svdir filesep 'macministat.bin'],'file') ~= 2
    
    disp('UPLOADCHECKER_PI: MacMinistat file not found. Creating now...');
    [fid, ~] = fopen([svdir filesep 'macministat.bin'],'w');
    fwrite(fid,timestamp,'double');
    fclose(fid);
    
else
    
    if vrb == 1, disp('UPLOADCHECKER_PI: Found macministat.bin! Writing timestamp.'); end
    try
        [mfid, ~] = fopen([svdir filesep 'macministat.bin'],'a');
        fwrite(mfid,timestamp,'double');
        fclose(mfid);
    catch
        disp('error opening macministat file on first try. Trying again.');
        fclose('all');
        [mfid, ~] = fopen([svdir filesep 'macministat.bin'],'a');
        fwrite(mfid,timestamp,'double');
        fclose(mfid);
    end
    
end




% TEST FOR READING THE FILE
%         fname = [];
%         fname = ['/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus/lastsavetime.bin'];
%         [fid, message] = fopen(fname,'r');
%         A = [];
%         A = fread(fid,'double');
%         fclose(fid);

