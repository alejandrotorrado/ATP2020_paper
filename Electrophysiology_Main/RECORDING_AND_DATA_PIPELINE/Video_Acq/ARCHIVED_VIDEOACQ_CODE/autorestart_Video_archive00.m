function autorestart_Video
% v2: ATP 11/23/15
% updated for new movement tracking
% added gui elements to choose framerate, and save frequency
%video timer GUI
imaqreset;
delete(timerfind);

global handles;

handles.XY_list_A=[];
handles.XY_list_B=[];
handles.XY_list=[];
handles.imcount = 0;
handles.viewIm = 1;

handles.piDir = '/Users/keithhengen/Desktop/RPI';
handles.checkerDir = [handles.piDir filesep 'RECORD'];
if ~exist(handles.checkerDir,'dir')
    disp('Rpi not found. Trying to mount it.');
    mountPi = system(['mount_smbfs //pi:invivo99@129.64.111.167/ExpLogs ' handles.piDir]);
    if mountPi ~= 0
        error('Could not mount RPi. Please check network connections.');
    else
        disp(['Successfully mounted Rpi in ' handles.piDir])
    end
else
    disp(['Found Rpi at ' handles.piDir])
end

disp('Loading video metadata.')

fid = fopen([handles.checkerDir filesep 'vid_metadata.txt']);
data = textscan(fid,'%s%s%s%s%s', 'HeaderLines', 1, 'CollectOutput', 0);
fclose(fid);

handles.reloaddata.animal = data{1}{1};
handles.reloaddata.date = data{2}{1};
savedir_anim = data{3}{1};
savedir = ['/Users/keithhengen/Desktop/animal video frames' filesep savedir_anim];
handles.reloaddata.savedir = savedir;
handles.reloaddata.framerate = str2double(data{4}{1});
handles.reloaddata.interval = str2double(data{5}{1});


disp(['Metadata is loaded. Working with animal ' handles.reloaddata.animal '.'])

disp('Initializing GUI.');

handles.vidf = figure('NumberTitle','off','Color','white','Visible','off',...
    'Position',[300 200 500 400],'DockControls','off','MenuBar','none',...
    'ToolBar','none','WindowStyle','modal');

handles.stopvid=uicontrol('Style','pushbutton','String','Stop Recording','Position',...
    [175 250 300 55],'FontSize',16,'BackgroundColor','red','CallBack',@stopvidacq);

handles.status_msg = uicontrol('style','text','string','Video acquisition restarted.',...
    'visible','off','Position',[175 150 300 80],'FontSize',16,'BackgroundColor','white');

set(handles.vidf,'visible','on');
set(handles.stopvid,'visible','on');


restartvidacq;



    function restartvidacq(~,~)
        handles.RUNSTATUS=1;
        set(handles.status_msg,'visible','on');
        % create mac mini restart flag
        vid_restart_time = unixtime(clock);
        vfid = fopen([handles.checkerDir filesep 'vid_restarted.bin'],'w');
        fwrite(vfid,vid_restart_time,'double');
        fclose(vfid);
        
        % delete mac mini crash flag
        if exist([handles.checkerDir filesep 'macminiCrash.bin'],'file') == 2
            delete([handles.checkerDir filesep 'macminiCrash.bin']);
        end
        
        video_run_step3_2animals(handles.reloaddata.animal,...
            handles.reloaddata.date,handles.reloaddata.savedir,...
            handles.reloaddata.framerate, handles.reloaddata.interval);
        figure(handles.vidf);
    end

    function stopvidacq(~,~)
        set(handles.stopvid,'Visible','off');
        set(handles.status_msg,'visible','off');
        handles.endvid=uicontrol('Style','text','BackgroundColor','green','position',...
            [175 250 400 55],'String','Video Finished. Close this window!',...
            'FontSize',16);
        [status,result] = system('python /Users/keithhengen/Desktop/killPyChecker.py');
        disp('status is:')
        disp(status)
        disp('result is:')
        disp(result)
        stop(handles.t);
        delete(handles.t);
        stop (handles.vid);
        delete (handles.vid);
        fclose('all');
        clear all;
        clc;

    end


end




