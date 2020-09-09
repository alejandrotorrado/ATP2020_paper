function VidTimerGui_2animals_v2
% v2: ATP 11/23/15
% updated for new movement tracking
% added gui elements to choose framerate, and save frequency
%video timer GUI
imaqreset;
delete(timerfind);

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

handles.XY_list_A=[];
handles.XY_list_B=[];
handles.XY_list=[];

handles.checkerDir = uigetdir(cd,'Pick the MacMini checker dropbox folder'); % temporary


f=figure('Visible','off','Color','white','Position',[300 200 800 700]);
ttext=uicontrol('Style','text','BackgroundColor','white','position',...
    [175 630 2*ctl_w ctl_h],'String','Set up time-stamped video of a recording.',...
    'FontSize',16);
A_num_text = uicontrol('Style','Text','BackgroundColor','white','position',...
    [175 570 ctl_w ctl_h],'String','How many animals?',...
    'FontSize',16);
AN_A_text=uicontrol('Style','text','BackgroundColor','white','position',...
    [A_left num_bott ctl_w ctl_h],'String','Animal A (left) Number:',...
    'FontSize',16);
AN_B_text=uicontrol('Style','text','BackgroundColor','white','position',...
    [B_left num_bott ctl_w ctl_h],'String','Animal B (right) Number:',...
    'FontSize',16);
AN_A_uitext=uicontrol('Style','edit','Position',[A_left num_bott-30 ctl_w ctl_h],...
    'CallBack',@animalAFunc);
AN_B_uitext=uicontrol('Style','edit','Position',[B_left num_bott-30 ctl_w ctl_h],...
    'CallBack',@animalBFunc);
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
A1button=uicontrol('Style','pushbutton','String','1','Position',...
    [420 585 ctl_h ctl_h],'FontSize',16,'CallBack',@Anum_1);
A2button=uicontrol('Style','pushbutton','String','2','Position',...
    [490 585 ctl_h ctl_h],'FontSize',16,'CallBack',@Anum_2);
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
        stopbutton=uicontrol('Style','pushbutton','String','Stop Recording','Position',...
            [550 125 ctl_w ctl_h],'FontSize',16,'BackgroundColor','red','CallBack',@stoprec_CallBack);
        handles.RUNSTATUS=1;
        
        
        % - - - - - - - - - 
        %    clear out the server log from this computer (if it's still
        %    there from an old recording)
        if exist('/Volumes/turrigiano-lab/RECORD/mministat.bin','file') == 2;
            delete('/Volumes/turrigiano-lab/RECORD/mministat.bin');
        end
        % - - - - - - - - - 
        
        if num_animals == 1
            a_num=get(AN_A_uitext,'String');
            a_date=get(D_uitext,'String');
            DR=handles.varA;
            video_run_step3(a_num,a_date,DR); % need to make this function
            figure (1);
        elseif num_animals == 2
            a_numA = get(AN_A_uitext,'String');
            a_numB = get(AN_B_uitext,'String');
            a_date=get(D_uitext,'String');
            DR_A = handles.varA;
%             DR_B = handles.varB;
            video_run_step3_2animals(a_numA,a_numB,a_date,DR_A);
            figure(1);
        else
            warning('Please select number of animals before running video');
        end
        
        function stoprec_CallBack(source,eventdata)
            set(stopbutton,'Visible','off');
            finish_text=uicontrol('Style','text','BackgroundColor','green','position',...
                [550 125 ctl_w ctl_h],'String','Video Finished. Close this window!',...
                'FontSize',16);
            stop(handles.t);
            delete(handles.t);
            stop (handles.vid);
            delete (handles.vid);
            delete(handles.ratTracker3000_A);
            delete(handles.ratTracker3000_B);
            clear all
            clc
        end
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
        

%{
    function calibratecamera(source, eventdata) 
        % Interact with a single frame from the camera's current position
        % and select the four corners of the image that will be searched
        % for an LED source.
        a = imaqhwinfo;
        [camera_name, camera_id, format] = getCameraInfo(a);
        vid = videoinput(camera_name, camera_id, format);
        triggerconfig(vid,'Manual')
        set(vid,'TriggerRepeat',0)
        set(vid,'FramesPerTrigger',1)
        
        start(vid)
        
        trigger(vid)
        data = getdata(vid,1);
        
        calIm = figure (6);
        set(gcf, 'MenuBar', 'none');
        set(gcf, 'ToolBar', 'none');
        set(gca,'XTick',[],'YTick',[]);
        imshow(data)
        
        %%%******%%%******************************************************%%%******%%%
        %%%******%%%___ put an IF statement using variable num_animals ___%%%******%%%
        %%%******%%%******************************************************%%%******%%%
        
        if num_animals == 1
            
            t1=text(650, 500, 'Click the top left limit','color','r','fontsize',24);
            [x1,y1]=ginput(1);
            delete(t1)
            t2=text(650, 500, 'Click the top right limit','color','r','fontsize',24);
            [x2,y2]=ginput(1);
            delete(t2)
            t3=text(650, 500, 'Click the bottom left limit','color','r','fontsize',24);
            [x3,y3]=ginput(1);
            delete(t3)
            t4=text(650, 500, 'Click the bottom right limit','color','r','fontsize',24);
            [x4,y4]=ginput(1);
            delete(calIm);
            
            handles.xStart = min([x1 x3]);
            handles.xEnd = max([x2 x4]);
            handles.yStart = min([y1 y2]);
            handles.yEnd = max([y3 y4]);
            
        elseif num_animals == 2
            
            t_A = text(500, 400, 'Frame on the left (Animal A)','color',[1.0 0.8 0.0],'fontsize',24);
            t1=text(500, 500, 'Click the top left limit','color','r','fontsize',24);
            [x1,y1]=ginput(1);
            delete(t1)
            t2=text(500, 500, 'Click the top right limit','color','r','fontsize',24);
            [x2,y2]=ginput(1);
            delete(t2)
            t3=text(500, 500, 'Click the bottom left limit','color','r','fontsize',24);
            [x3,y3]=ginput(1);
            delete(t3)
            t4=text(500, 500, 'Click the bottom right limit','color','r','fontsize',24);
            [x4,y4]=ginput(1);
            delete(t4);
            delete(t_A);
            
            t_B = text(1300, 400, 'Frame on the right (animal B)','color',[1.0 0.8 0.0],'fontsize',24);
            t1=text(1300, 500, 'Click the top left limit','color','r','fontsize',24);
            [x5,y5]=ginput(1);
            delete(t1)
            t2=text(1300, 500, 'Click the top right limit','color','r','fontsize',24);
            [x6,y6]=ginput(1);
            delete(t2)
            t3=text(1300, 500, 'Click the bottom left limit','color','r','fontsize',24);
            [x7,y7]=ginput(1);
            delete(t3)
            t4=text(1300, 500, 'Click the bottom right limit','color','r','fontsize',24);
            [x8,y8]=ginput(1);
            delete(t4);
            delete(t_B);
            
            delete(calIm);
            
            handles.xStart_A = min([x1 x3]);
            handles.xEnd_A = max([x2 x4]);
            handles.yStart_A = min([y1 y2]);
            handles.yEnd_A = max([y3 y4]);
            
            handles.xStart_B = min([x5 x7]);
            handles.xEnd_B = max([x6 x8]);
            handles.yStart_B = min([y5 y6]);
            handles.yEnd_B = max([y7 y8]);
            
        else
            imwarn=warndlg('Select number of animals first.',...
                'Do it right, son.','modal');
            uiwait(imwarn);
            delete(calIm);
            delete(f);
        end
        
    end
%}
end




