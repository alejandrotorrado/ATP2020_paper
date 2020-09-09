function uploadchecker(T_seconds)

% UPLOADCHECKER: Create and save a binary file containing the unix time
% stamp of the last image recorded by the rat video recording software.
% Save this data to drop box such that another computer linked to the
% dropbox account can check to see if the saved variable is increasing
% every 30s or so. If it is not, then the user can deduce that the mac mini
% (responsible for recording video) has crashed. This can then trigger a
% text message to be sent to the user.

if ~exist('/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus','dir');
    mkdir('/Users/keithhengen/Dropbox/VIDEOTRACK','macministatus');
end


if ~exist('/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus/lastsavetime.bin','file');
  
    [fid, message] = fopen('/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus/lastsavetime.bin','a');
    fwrite(fid,T_seconds,'double')
    fclose(fid);
    
else
    
    [fid, message] = fopen('/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus/lastsavetime.bin','w');
    fwrite(fid,T_seconds,'double');
    fclose(fid);
    
end




% TEST FOR READING THE FILE
%         fname = [];
%         fname = ['/Users/keithhengen/Dropbox/VIDEOTRACK/macministatus/lastsavetime.bin'];
%         [fid, message] = fopen(fname,'r');
%         A = [];
%         A = fread(fid,'double');
%         fclose(fid);

