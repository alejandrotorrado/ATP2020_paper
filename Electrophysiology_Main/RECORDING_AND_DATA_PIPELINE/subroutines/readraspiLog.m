function [onTimes offTimes C] =  readraspiLog(fname,recgo)
% This function will load the file indicated by 'fname', which should
% correspond to a log file created during one of our (KBH, ATP) recordings
% and saved on the raspberry pi server. The log file contains data in the
% format of rows such as:
%
%        '2016-03-05 09:05:10	Started KH6364 experiment.'
% 
% The following code will read the timestamp of each of the log events as
% well as the string statement associated with it. Additionally, this will
% recognize key events as starts and stops in the data acquisition, and
% will return a matrix indicating those times. This matrix should be
% applied to all cells from relevant animals to have highly accurate
% information about when cells are on/off line (not considering the UI
% sensitive determination of the cell's presence and health).


fid = fopen(fname);

A = textscan(fid,'%u %u %u %*[^\n]','Delimiter','-');

fseek(fid,0,'bof');
B = textscan(fid,'%*[^ ] %u %u %u %*[^\n]','delimiter',':');

fseek(fid,0,'bof');
C = textscan(fid,'%*[^\t] %s %*[^\n]','delimiter','\t');
C = C{1};

fclose(fid);

dts = double([A{1} A{2} A{3} B{1} B{2} B{3}]);
udts = unixtime(dts);


mtx(:,1) = udts;
mtx(:,2) = zeros;
% Index the rows that contain specific language:
% start/stop experiment
expGO           =  ~cellfun( @isempty,( regexp(C,'Started(\s*\w*\s*[Experiment])') ));
mtx(expGO,2)    = 101;
MLSTOP          =  ~cellfun( @isempty,( regexp(C,'MATLAB: Stopped recording for(s*\w*)')));
mtx(MLSTOP,2)   = 102;

% software dies:
SOFTCRASH       =  ~cellfun( @isempty,( regexp(C,'Found software crash in loopchecker. Created softwareCrash flag.')));
mtx(SOFTCRASH,2)= 201;
deadSOFT        =  ~cellfun( @isempty,( regexp(C,'Executing startover. Killing matlab and synapse.')));
mtx(deadSOFT,2) = 202;
rstSOFT         =  ~cellfun( @isempty,( regexp(C,'Executing startover. Restarting matlab and synapse.')));
mtx(rstSOFT,2)  = 203;

% RS4 crash:
RS4CRASH        =  ~cellfun( @isempty,( regexp(C,'Found RS4 crash in loopchecker. Creating rs4flag.')));
mtx(RS4CRASH,2) = 301;
RS4SOFT         =  ~cellfun( @isempty,( regexp(C,'Killing software following RS4 crash. Killing matlab and synapse.')));
mtx(RS4SOFT,2)  = 302;

% general restart
% rstREC       =  ~cellfun( @isempty,( regexp(C,'Found restart flag in startover. Recording restarted.' ) ) ) ;
% mtx(rstREC,2)= 401;
reboot       =  ~cellfun( @isempty,( regexp(C,'Could not find restart flag after 60 sec. Rebooting startover.')));
mtx(reboot,2)= 402;
flgGO        =  ~cellfun( @isempty,( regexp(C,'Found restart flag in startover. Recording restarted.')));
mtx(flgGO,2) = 403;

% software up and running
PYGO          =  ~cellfun( @isempty,( regexp(C,'Starting main python rig code to monitor recording.')));
mtx(PYGO,2)   = 501;
LOOPGO        =  ~cellfun( @isempty,( regexp(C,'Found rigstat flag in main python rig code. Entering loopchecker.')));
mtx(LOOPGO,2) = 502;
LOOPit        =  ~cellfun( @isempty,( regexp(C,'Resuming loopchecker.')));
mtx(LOOPit,2) = 503;

% stop recording
ENDREC        =  ~cellfun( @isempty,( regexp(C,'Found STOPREC flag in loopchecker. Quitting python.')));
mtx(ENDREC,2) = 601;

% Find Recording Starts:
start_seq = [101 501 502]; 
startIDX  = strfind(mtx(:,2)',start_seq);

st__over_seq = [202 203 403 503]; 
st_overIDX   = strfind(mtx(:,2)',st__over_seq) + 2;

sft_crash_seq = 201;
sft_crashIDX  = strfind(mtx(:,2)',sft_crash_seq);

RS4_crash_seq = [301 302];
RS4_crashIDX  = strfind(mtx(:,2)',RS4_crash_seq);

endrec_seq = 601;
endrecIDX  = strfind(mtx(:,2)',endrec_seq);

mtl_stop_seq = 102;
mtl_stopIDX  = strfind(mtx(:,2)',mtl_stop_seq);

onT(:,1) = mtx([startIDX st_overIDX], 1);
onT(:,2) = ones;

offT(:,1) = mtx([sft_crashIDX RS4_crashIDX endrecIDX mtl_stopIDX],1);
offT(:,2) = zeros;

onOff = sortrows([onT; offT],1);

% Add line that if last two entries in onOff are zeros that the second one
% can be ignored. That's the rig shutting down, effectively.
if strfind(onOff(end-1:end,2)', [0 0]);
    onOff(end,:) = [];
end
    
% adjust start times to align on 7:30 am on the first day of the recording
tdiff = recgo  - 27000;
trem = rem(tdiff,(24*3600));

onOff(:,1) = onOff(:,1) - recgo + trem;

% Separate on and off times now for return to caller.
onTimes     = onOff(onOff(:,2)==1,:);
offTimes    = onOff(onOff(:,2)==0,:);











