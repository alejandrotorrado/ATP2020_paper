function VidTimerGui_2animals_v3
% v2: ATP 11/23/15
% updated for new movement tracking
% added gui elements to choose framerate, and save frequency
%video timer GUI
imaqreset;
delete(timerfind);
fclose('all');
clear all;
clc;



% define standard position and size for GUI controls
ctl_h = 55; % height
ctl_w = 200; % width
A_left = 150;
B_left = 400;
num_bott = 480;
dir_bott = 360;
date_bott = 280;
num_animals = 0;

global handles;

% set to 0 to suppress output
handles.imcount = 0;
handles.viewIm = 1;
if handles.viewIm
    disp('Current settings will display camera output.')
    noOut = input('Do you want to suppress camera output?  (y/n):  ','s');
    if strcmp(noOut,'y')
        handles.viewIm = 0;
        disp('Camera output suppressed.');
    elseif strcmp(noOut,'n')
        handles.viewIm = 1;
        disp('Camera ouptut will be displayed.');
    else
        error('Invalid input.');
    end
end

handles.XY_list_A=[];
handles.XY_list_B=[];
handles.XY_list=[];

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


f=figure('Visible','off','Color','white','Position',[300 200 800 700]);
ttext=uicontrol('Style','text','BackgroundColor','white','position',...
    [175 630 2*ctl_w ctl_h],'String','Set up time-stamped video of a recording.',...
    'FontSize',16);
% A_num_text = uicontrol('Style','Text','BackgroundColor','white','position',...
%     [175 570 ctl_w ctl_h],'String','How many animals?',...
%     'FontSize',16);
AN_A_text=uicontrol('Style','text','BackgroundColor','white','position',...
    [175 num_bott 2*ctl_w ctl_h],'String','Animal name:',...
    'FontSize',16);
% AN_B_text=uicontrol('Style','text','BackgroundColor','white','position',...
%     [B_left num_bott ctl_w ctl_h],'String','Animal B (right) Number:',...
%     'FontSize',16);
AN_A_uitext=uicontrol('Style','edit','Position',[175 num_bott-30 2*ctl_w ctl_h],...
    'CallBack',@animalAFunc);
% AN_B_uitext=uicontrol('Style','edit','Position',[B_left num_bott-30 ctl_w ctl_h],...
%     'CallBack',@animalBFunc);
D_text=uicontrol('Style','text','BackgroundColor','white','position',...
    [175 date_bott 2*ctl_w ctl_h],'String','Date (mmddyy):',...
    'FontSize',16);
D_uitext=uicontrol('Style','edit','Position',[175 date_bott-30 2*ctl_w ctl_h],...
    'CallBack',@dateFunc);
FR_text = uicontrol('Style','text','BackgroundColor','white','position',...
    [A_left date_bott-110 ctl_w ctl_h],'String','Framerate (Hz):',...
    'FontSize',16);
Int_text = uicontrol('Style','text','BackgroundColor','white','position',...
    [B_left date_bott-110 ctl_w ctl_h],'String','Save interval (sec):',...
    'FontSize',16);
FR_uitext = uicontrol('Style','edit','position',...
    [A_left date_bott-140 ctl_w ctl_h],'Callback',@frameFunc);
Int_uitext = uicontrol('Style','edit','position',...
    [B_left date_bott-140 ctl_w ctl_h],'Callback',@intervalFunc);

handles.varA='empty';
handles.varB='empty';

%push button to call uigetdir - find a directory to save the video
%acquisition
% A1button=uicontrol('Style','pushbutton','String','1','Position',...
%     [420 585 ctl_h ctl_h],'FontSize',16,'CallBack',@Anum_1);
% A2button=uicontrol('Style','pushbutton','String','2','Position',...
%     [490 585 ctl_h ctl_h],'FontSize',16,'CallBack',@Anum_2);
DAbutton=uicontrol('Style','pushbutton','String','Save Directory','Position',...
    [225 dir_bott 1.5*ctl_w ctl_h],'FontSize',16,'CallBack',@pickdirA);
% DBbutton=uicontrol('Style','pushbutton','String','Animal B Directory','Position',...
%     [B_left dir_bott ctl_w ctl_h],'FontSize',16,'CallBack',@pickdirB);
Rbutton=uicontrol('Style','pushbutton','String','Start Recording','Position',...
    [550 70 ctl_w ctl_h],'FontSize',16,'BackgroundColor','red','CallBack',@runvidacq);
PULSEbutton=uicontrol('Style','pushbutton','String','Test Pulses?','Position',...
    [275 70 ctl_w ctl_h],'FontSize',16,'CallBack',@pulsetest,'visible','on');
% CALIBRATEbutton=uicontrol('Style','pushbutton','String','Calibrate Camera','Position',...
%     [300 90 ctl_w ctl_h],'FontSize',16,'CallBack',@calibratecamera,'visible','on');


%Now the GUI is made visible
set(f,'Visible','on');


%animal number text entry function
    function animalAFunc(source,eventdata)
        %called by pressing enter in the animal number text field
        set(AN_A_uitext,'Visible','off');
        print_animal=get(AN_A_uitext,'String');
        show_animal=uicontrol('Style','text','BackgroundColor','white',...
            'Position',[A_left num_bott-40 ctl_w ctl_h],'String',print_animal,'FontSize',...
            18,'ForegroundColor','red');
        set(show_animal,'Visible','on')
    end

    function animalBFunc(source,eventdata)
        %called by pressing enter in the animal number text field
        set(AN_B_uitext,'Visible','off');
        print_animal=get(AN_B_uitext,'String');
        show_animal=uicontrol('Style','text','BackgroundColor','white',...
            'Position',[B_left num_bott-40 ctl_w ctl_h],'String',print_animal,'FontSize',...
            18,'ForegroundColor','red');
        set(show_animal,'Visible','on')
    end

% choose number of animals callback function
    function Anum_1(source,eventdata)
        set([AN_B_text, AN_B_uitext,DBbutton],'Visible','off');
        set([A1button, A2button],'Visible','off');
        num_animals = 1;
        print_an_num = 'Using 1 animal';
        show_an_num = uicontrol('Style','Text','BackgroundColor','white',...
            'Position',[420,570,ctl_w,ctl_h],'String',print_an_num,...
            'FontSize',18,'ForegroundColor','red');
        set(show_an_num,'Visible','on');
    end

    function Anum_2(source,eventdata)
        set([A1button, A2button],'Visible','off');
        print_an_num = 'Using 2 animals';
        num_animals = 2;
        show_an_num = uicontrol('Style','Text','BackgroundColor','white',...
            'Position',[420,570,ctl_w,ctl_h],'String',print_an_num,...
            'FontSize',18,'ForegroundColor','red');
        set(show_an_num,'Visible','on');
    end

%date text entry function
    function dateFunc(source,eventdata)
        %called by pressing enter in the date text field
        set(D_uitext,'Visible','off');
        print_date=get(D_uitext,'String');
        show_date=uicontrol('Style','text','BackgroundColor','white',...
            'Position',[175 date_bott-30 400 ctl_h],'String',print_date,'FontSize',...
            18,'ForegroundColor','red');
        set(show_date,'Visible','on')
    end

%Directory button function
    function pickdirA(source,eventdata)
        %called by Dbutton press
        set(DAbutton,'Visible','off');
        save_dir=uigetdir;
        handles.varA=save_dir;
        dir_str=uicontrol('Style','text','BackgroundColor','white','Position',...
            [200 dir_bott-20 1.75*ctl_w ctl_h],'String',save_dir,'ForegroundColor','red',...
            'FontSize',12,'HorizontalAlignment','left');
        stxt=uicontrol('Style','text','BackgroundColor','white','position',...
            [A_left-50 dir_bott+40 0.75*ctl_w 0.5*ctl_h],'String','Saving to:',...
            'FontSize',12,'Visible','off','ForegroundColor','black',...
            'HorizontalAlignment','center');
        set([stxt dir_str],'Visible','on')
        
    end

%     function pickdirB(source,eventdata)
%         %called by Dbutton press
%         set(DBbutton,'Visible','off');
%         save_dir=uigetdir;
%         handles.varB=save_dir;
%         dir_str=uicontrol('Style','text','BackgroundColor','white','Position',...
%             [B_left+50 dir_bott-20 1.5*ctl_w ctl_h],'String',save_dir,'ForegroundColor','red',...
%             'FontSize',10,'HorizontalAlignment','left');
%         stxt=uicontrol('Style','text','BackgroundColor','white','position',...
%             [B_left+50 dir_bott+40 0.75*ctl_w 0.5*ctl_h],'String','Saving to:',...
%             'FontSize',10,'Visible','off','ForegroundColor','black',...
%             'HorizontalAlignment','left');
%         set([stxt dir_str],'Visible','on')
%
%     end

    function runvidacq(source,eventdata) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change here for 2 animals
        set(Rbutton,'Visible','off');
        handles.stopbutton=uicontrol('Style','pushbutton','String','Stop Recording','Position',...
            [550 125 ctl_w ctl_h],'FontSize',16,'BackgroundColor','red','CallBack',@stoprec_CallBack);
        handles.RUNSTATUS=1;
        
        
        % - - - - - - - - -
        %    clear out the macministatus file and video metadata (if it's still
        %    there from an old recording)
        if exist([handles.checkerDir filesep 'macministat.bin'],'file') == 2;
            delete([handles.checkerDir filesep 'macministat.bin']);
        end
        if exist([handles.checkerDir filesep 'vid_metadata.txt'],'file') == 2;
            delete([handles.checkerDir filesep 'vid_metadata.txt']);
        end
        % - - - - - - - - -
        a_num=get(AN_A_uitext,'String');
        a_date=get(D_uitext,'String');
        
        % python monitor file name
        pyfname = ['checkmacmini_' a_num '.txt'];
        % run the macmini python code that will monitor macministat
        pycmd = ['python /Users/keithhengen/Desktop/checkMacmini.py > /Users/keithhengen/Desktop/' pyfname ' &'];
        [status,result] = system(pycmd);
        
        %         if num_animals == 1
       
        %             DR=handles.varA;
        %             video_run_step3(a_num,a_date,DR); % need to make this function
        %             figure (1);
        %         elseif num_animals == 2
        %             a_numA = get(AN_A_uitext,'String');
        %             a_numB = get(AN_B_uitext,'String');
        %             a_date=get(D_uitext,'String');
        DR_A = handles.varA;
        %             DR_B = handles.varB;
        video_run_step3_2animals(a_num,a_date,DR_A,handles.framerate,handles.interval);
        figure(f);
        %         else
        %             warning('Please select number of animals before running video');
    end

    function stoprec_CallBack(source,eventdata)
        set(handles.stopbutton,'Visible','off');
        finish_text=uicontrol('Style','text','BackgroundColor','green','position',...
            [550 125 ctl_w ctl_h],'String','Video Finished. Close this window!',...
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


    function pulsetest (source, eventdata)
        
        StimTriggerOpen;
        
        ii=0;
        while ii<=50;
            ii=ii+1;
            
            StimTriggerAct('Pin35On');
            pause(0.25)
            StimTriggerAct('Pin35Off');
            pause(0.25);
            disp(['Pulse ' num2str(ii) ' sent.']);
        end
        
    end

    function frameFunc (source,eventdata)
        %called by pressing enter in the date text field
        set(FR_uitext,'Visible','off');
        print_FR=get(FR_uitext,'String');
        show_FR=uicontrol('Style','text','BackgroundColor','white',...
            'Position',[A_left date_bott-140 ctl_w ctl_h],'String',...
            [print_FR ' Hz'],'FontSize',16,'ForegroundColor','red');
        set(show_FR,'Visible','on')
        handles.framerate = str2double(print_FR);
    end

    function intervalFunc (source,eventdata)
        %called by pressing enter in the date text field
        set(Int_uitext,'Visible','off');
        print_Int=get(Int_uitext,'String');
        show_Int=uicontrol('Style','text','BackgroundColor','white',...
            'Position',[B_left date_bott-140 ctl_w ctl_h],'String',...
            ['Save every ' print_Int ' seconds'],'FontSize',16,'ForegroundColor','red');
        set(show_Int,'Visible','on')
        handles.interval = str2double(print_Int);
    end

end






