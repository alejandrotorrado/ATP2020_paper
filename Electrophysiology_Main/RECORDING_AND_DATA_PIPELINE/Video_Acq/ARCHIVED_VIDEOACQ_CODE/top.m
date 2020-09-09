
% NEW VIDEO CODE - COMPLETELY REWRITTEN ON 05.07.16 KBH. NEEDS WORK. BE
% CAREFUL...

clear all
close all
clc
global handles

a_num = 'dryrun';
today = '052316';
saveinterval = 1800;

handles.vFile = [a_num '.avi'];
handles.vDir = ['/Users/keithhengen/Desktop/animal_video_frames' filesep a_num];

if ~exist(handles.vDir,'dir');
    mkdir(handles.vDir);
end

handles.v = VideoWriter([handles.vDir handles.vFile]);

handles.interval = 0.5;
handles.fileCount = 0;
handles.frametimes = [];
handles.cam  = webcam;
handles.t = timer('StartDelay', 0, 'Period', handles.interval, 'TasksToExecute', inf, ...
    'ExecutionMode', 'fixedRate');
handles.t.StartFcn = {@AnnounceStartVid, 'Starting video timer'};
handles.t.TimerFcn = {@bottom, 'My message', saveinterval};
handles.count = 0;

% - - - - Check to see if there are already video files in place. If so,
% find the highest number and then start adding to that. Otherwise, you'll
% overwrite your data.
backcheck = dir([handles.vDir filesep '*.avi']);
backnames = {backcheck.name};

pattern = '\d+[_]+\d*+[_]'; ncount = 0;
for ee = 1:size(backnames,2);
    test = regexp(backnames{ee},pattern);
    if ~isempty(test);
        ncount = ncount + 1;
        t2 = regexp(backnames{ee},'_');
        nameno(ncount) = str2double(backnames{ee}(t2(2)+1:end-4));
    end
end
if exist('nameno','var');
    handles.fileCount = max(nameno);
else
    handles.fileCount = 0;
end
% - - - - - - - - - - - - - - - - - - 



% - - - - - deal with raspi - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - 
% - - - - - - - - - - - - mount the pi - - - - - - - - - - - - - - - - - - 
% handles.piDir = '/Users/keithhengen/Desktop/RPI';
% handles.checkerDir = [handles.piDir filesep 'RECORD'];
% if ~exist(handles.checkerDir,'dir')
%     disp('Rpi not found. Trying to mount it.');
%     mountPi = system(['mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs ' handles.piDir]);
%     if mountPi ~= 0
%         error('Could not mount RPi. Please check network connections.');
%     else
%         disp(['Successfully mounted Rpi in ' handles.piDir])
%     end
% else
%     disp(['Found Rpi at ' handles.piDir])
% end
% % - - - - - - - 
% 
% if exist([handles.checkerDir filesep 'macministat.bin'],'file') == 2;
%     delete([handles.checkerDir filesep 'macministat.bin']);
% end
% if exist([handles.checkerDir filesep 'vid_metadata.txt'],'file') == 2;
%     delete([handles.checkerDir filesep 'vid_metadata.txt']);
% end
% - - - - - - - - -
%a_num=get(AN_A_uitext,'String');
%a_date=get(D_uitext,'String');

% % python monitor file name
% pyfname = ['checkmacmini_' a_num '.txt'];
% % run the macmini python code that will monitor macministat
% pycmd = ['python /Users/keithhengen/Desktop/checkMacmini.py > /Users/keithhengen/Desktop/' pyfname ' &'];
% [status,result] = system(pycmd);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Do the metadata
direct = handles.vDir;
handles.timezero = unixtime(clock);
handles.FILESAVEDIRECTORY=direct;
slashes = regexp(direct,filesep);
handles.FILESAVEDIR_ANIM = a_num; %direct(slashes(end)+1:end);

handles.ANIMAL_NUMBER=a_num;
% handles.ANIMAL_B_NUMBER=animalB;
handles.TODAYS_DATE=today;

% save video metadata file if not there - - - this simply holds the
% following information: animal, date, save direc, frame rate, save
% interval. The top line (first frprint) writes the header titles. The next
% fprint line writes the data in place.
% if exist([handles.checkerDir filesep 'vid_metadata.txt'],'file') ~= 2
%     fid = fopen([handles.checkerDir filesep 'vid_metadata.txt'],'w');
%     fprintf(fid,'%s\t %s\t %s\t %s\t %s\n','animal','date','savedir','framerate','save interval');
%     fprintf(fid,'%s\t %s\t %s\t %s\t %s\n',a_num,today,handles.FILESAVEDIR_ANIM,num2str(1/handles.interval),saveinterval);
%     fclose(fid);
% else
%     disp('Found video metadata file.')
% end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - -

% - - - - - 

open(handles.v)
start(handles.t)



% KILL THE PYTHON CODE:
%[status,result] = system('python /Users/keithhengen/Desktop/killPyChecker.py')




