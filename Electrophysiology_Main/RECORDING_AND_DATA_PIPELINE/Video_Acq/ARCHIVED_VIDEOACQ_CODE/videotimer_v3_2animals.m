function videotimer_v3_2animals(framerate,saveInterval)
% v3: ATP 11/23/2015
% Updating for new movement tracking.
%  1) removing everything related to LED tracking
%  2) Only need one saving directory, because same video is used for both
%     animals; ROI definition will come later, in movement tracking code
%  3) Increase frame rate
%  4) Keep timing pulses
%
% First pass: comment out old stuff that is not needed; check if this
% works; then delete and make code final version.
%
% old info:
%   for use with animal behavior video-collection scripts. this is the secondary
%   function in the line-up. this sets up a timer and collects global variables re: data
%   saving directory,animal number, and date. Finally, script starts the timer
%
% global TODAYS_DATE;
% global ANIMAL_NUMBER;
% global FILESAVEDIRECTORY;
% FILESAVEDIRECTORY=uigetdir('Pick a directory for saving the data');
% ANIMAL_NUMBER=input('What is the animal number?','s');
% TODAYS_DATE=input('What is the date today? (mmddyy)','s');

global handles;

handles.timezero = unixtime(clock);

handles.saveframe=[];
% handles.XY_list_A=[];
% handles.XY_list_B=[];
% handles.smallcount_A=0;
% handles.smallcount_B=0;
handles.framecount=0;
handles.resetTracker=0;

hDir_A = dir(handles.FILESAVEDIRECTORY);
% hDir_B = dir(handles.FILESAVEDIRECTORY_B);

cnt_A = 0; killem_A = [];
cnt_B = 0; killem_B = [];

for ii = 1:length(hDir_A);
    temp = [];
    temp = regexp(hDir_A(ii).name,'.mat', 'once');
    if isempty(temp);
        cnt_A = cnt_A+1;
        killem_A(cnt_A) = ii;
    end
end
hDir_A(killem_A) = [];
% 
% for ii = 1:length(hDir_B);
%     temp = [];
%     temp = regexp(hDir_B(ii).name,'.mat', 'once');
%     if isempty(temp);
%         cnt_B = cnt_B+1;
%         killem_B(cnt_B) = ii;
%     end
% end
% hDir_B(killem_B) = [];

if ~isempty(hDir_A);
    for ii = 1:length(hDir_A);
        temp1 = []; temp2 = [];
        temp1 = regexp(hDir_A(ii).name,'_');
        temp2 = regexp(hDir_A(ii).name,'.mat');
        fileList_A(ii,1) = str2double(  hDir_A(ii).name (temp1+1 : temp2-1)  );
        fileList_A(ii,2) = ii;
    end
end
% 
% if ~isempty(hDir_B);
%     for ii = 1:length(hDir_B);
%         temp1 = []; temp2 = [];
%         temp1 = regexp(hDir_B(ii).name,'_');
%         temp2 = regexp(hDir_B(ii).name,'.mat');
%         fileList_B(ii,1) = str2double(  hDir_B(ii).name (temp1+1 : temp2-1)  );
%         fileList_B(ii,2) = ii;
%     end
% end

if exist('fileList_A','var');
    handles.count_A = max(fileList_A(:,1)) + 1;
    disp(['detected other files in the save-to directory. count is starting at ' num2str(handles.count_A)]);
else
    handles.count_A = 0;
    disp(['did not detect other files in the save-to directory. count is starting at ' num2str(handles.count_A)]);
end
% 
% if exist('fileList_B','var');
%     handles.count_B = max(fileList_B(:,1)) + 1;
%     disp(['detected other files in the save-to directory. count is starting at ' num2str(handles.count_B)]);
% else
%     handles.count_B = 0;
%     disp(['did not detect other files in the save-to directory. count is starting at ' num2str(handles.count_B)]);
% end

% ATP 11/23/15
% In new code, 'Period' property in the timer object will determine the
% frame rate, as a frame is saved everytime the timer object is called. The
% period is the interval between frame saves.
% The framerate is set in the starting GUI and stored in handles. Calculate
% corresponding period (to 2 decimal places) and store in handles.
% handles.period = 1/handles.framerate;
disp(['frame rate ' num2str(double(framerate))]);
try
handles.t = timer('StartDelay', 0, 'Period', 1/double(framerate), 'TasksToExecute', inf, ...
    'ExecutionMode', 'fixedRate');
catch
    keyboard
end
disp(['period: ' num2str(1/double(framerate))]);
try
handles.t.StartFcn = {@my_callback_fcn, 'Starting video timer'};

handles.t.TimerFcn = {@timeVid_v4_2animals, 'My message', saveInterval};
catch
    keyboard
end


start(handles.t)
