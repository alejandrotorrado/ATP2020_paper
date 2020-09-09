function video_run_step3(animalN,today,direcT)
%This is the primary script for animal behavior video acquisition. The 
% principal goal is to have synchronized video acquisition and 
% electrophysiology data. This sets up the video input, starts the camera, 
% and then initializes the digitalinput/output through the parallel port. The
% parallel port is used to send a binary signal  to a data acquisition system
% every time a frame is collected. this binary signal is an "ID tag" that can
% be matched to the time of data collection on the acquitisition system.

%set up global variables that will be shared amongst the necessary scripts

global handles;

handles.FILESAVEDIRECTORY=direcT;
handles.ANIMAL_NUMBER=animalN;
handles.TODAYS_DATE=today;


a = imaqhwinfo;
[camera_name, camera_id, format] = getCameraInfo(a);
handles.vid = videoinput(camera_name, camera_id, format);
triggerconfig(handles.vid,'Manual')
set(handles.vid,'TriggerRepeat',inf)
set(handles.vid,'FramesPerTrigger',1)

start(handles.vid)

handles.count_A = 0;


StimTriggerOpen;

% THIS DOESN'T WORK ON A MAC. IS WRONG:
% dio = digitalio('parallel')
% addline(dio,7,0,'out')   % add pin 9 as line out
% global uddobj;
% uddobj = daqgetfield(dio,'uddobject') %uddobject is a variable within the dio structure. this is fast.

%call the secondary script
pause(1);
videotimer_v3_2animals
try
wait(handles.vid,1,'logging'); 
catch
%     disp('caught at 44 in video_run_step2.m');
%     keyboard
end

end