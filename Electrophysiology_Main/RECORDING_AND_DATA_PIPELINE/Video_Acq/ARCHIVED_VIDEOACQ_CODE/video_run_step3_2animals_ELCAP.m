function video_run_step3_2animals(animal,today,direct,frate,intrv)
%v3: ATP 11/23/2015
% Changing this for new movement tracking.
%
%This is the primary script for animal behavior video acquisition. The 
% principal goal is to have synchronized video acquisition and 
% electrophysiology data. This sets up the video input, starts the camera, 
% and then initializes the digitalinput/output through the parallel port. The
% parallel port is used to send a binary signal  to a data acquisition system
% every time a frame is collected. this binary signal is an "ID tag" that can
% be matched to the time of data collection on the acquitisition system.

%set up global variables that will be shared amongst the necessary scripts

global handles;
handles.timezero = unixtime(clock);
handles.FILESAVEDIRECTORY=direct;
slashes = regexp(direct,filesep);
handles.FILESAVEDIR_ANIM = direct(slashes(end)+1:end);

handles.ANIMAL_NUMBER=animal;
% handles.ANIMAL_B_NUMBER=animalB;
handles.TODAYS_DATE=today;

% save video metadata file if not there
if exist([handles.checkerDir filesep 'vid_metadata.txt'],'file') ~= 2
    fid = fopen([handles.checkerDir filesep 'vid_metadata.txt'],'w');
    fprintf(fid,'%s\t %s\t %s\t %s\t %s\n','animal','date','savedir','framerate','save interval');
    fprintf(fid,'%s\t %s\t %s\t %s\t %s\n',animal,today,handles.FILESAVEDIR_ANIM,num2str(frate),num2str(intrv));
    fclose(fid);
else
    disp('Found video metadata file.')
end
keyboard
a = imaqhwinfo;
[camera_name, camera_id, format] = getCameraInfo(a);
handles.vid = videoinput(camera_name, camera_id, format);
handles.vid = videoinput(b{1}, camera_name);

triggerconfig(handles.vid,'Manual')
set(handles.vid,'TriggerRepeat',inf)
set(handles.vid,'FramesPerTrigger',1)
set(handles.vid,'Timeout',50);

start(handles.vid)

handles.count_master=0;
handles.count_A=0;
% handles.count_B=0;



StimTriggerOpen;

% THIS DOESN'T WORK ON A MAC. IS WRONG:
% dio = digitalio('parallel')
% addline(dio,7,0,'out')   % add pin 9 as line out
% global uddobj;
% uddobj = daqgetfield(dio,'uddobject') %uddobject is a variable within the dio structure. this is fast.

%call the secondary script
pause(1);
videotimer_v3_2animals(frate,intrv)
try
wait(handles.vid,1,'logging'); 
catch
%     disp('caught at 44 in video_run_step2.m');
%     keyboard
end

end