%% BEHAVIORAL CODING GUI
%
% POSSIBLE IMPROVEMENT TO THIS:
% Random Forest (RF) may not be the best choice for this problem.
% Maybe look into Error Correcting Output Codes (ECOC)

function AutoSleepWakeCode_GUI_MLS_PyBETA_QuietEdit
warning('off');
global myname
startup
% set environment variable to include anaconda's python libraries.
% this will depend on the machine used. The global variable myname (set in
% the startup script) is used to figure this out.
if ismac
    pyPATH = add_anaconda_to_path();
    disp(pyPATH);
elseif ispc
    if strcmp(myname,'Jerboa')
        %         addpy = 'C:\Python27';
        %         addpy = 'C:\Python37-32';
        addpy =  'C:\Users\laneb\AppData\Local\Programs\Python\Python37-32';
    elseif strcmp(myname,'marmoset')
        addpy = '';
    elseif strcmp(myname,'Malabar')
        addpy = 'C:\Program Files\Python36';
    else
        addpy = '';
    end
    pyPATH = add_anaconda_to_path(addpy);
    disp(pyPATH);
end
setenv('PATH',pyPATH);

disp('Initializing GUI');

set(0,'Units','normalized')

%% GUI elements
% muscle beach logo
if ismac
    mb_logo = imread('/Volumes/turrigiano-lab/ATP_MAIN/misc_stuff/musclebeach_logo.tiff');
elseif ispc
    mb_logo = imread('Z:\ATP_MAIN\misc_stuff\musclebeach_logo.tiff');
end

% background color
f_color = [189, 243, 250]./255;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Main figure and loading buttons

% main window
handles.f = figure('WindowStyle','normal','HitTest','on','Visible','on','NumberTitle',...
    'off','Color',f_color,'Units','normalized','Position',[.05 .1 .9 .8]);

% "Load prior data?" along with "Yes" and "No" buttons
handles.priortext = uicontrol('Style','text','visible','on','Backgroundcolor',f_color,...
    'String','Load prior data?','units','normalized','Position',[0.42 0.83 0.1 0.04],...
    'Horizontalalignment','center','FontSize',16,'fontname','Times');

handles.loadprioryes = uicontrol('Style','pushbutton','String','Yes',...
    'units','normalized','Position', [0.41 0.78 0.05 0.04],'Callback',@LoadPrior,...
    'fontname','Times','FontSize',15,'Visible','on');

handles.loadpriorno = uicontrol('Style','pushbutton','String','No',...
    'units','normalized','Position', [0.49 0.78 0.05 0.04],'Callback',@LoadPrior,...
    'fontname','Times','FontSize',15,'Visible','on');

% "Load LFP" button
handles.loadLFP_button = uicontrol('Style','pushbutton','String','Load LFP data',...
    'units','normalized','Position', [0.42 0.58 0.2 0.1],'Callback',@LoadLFPdata,...
    'fontname','Times','FontSize',16,'Visible','off');

% "Have EMG?" "yes" and "no"
handles.useEMG_text = uicontrol('Style','text','String','Have EMG data?',...
    'units','normalized','Position', [0.42 0.58 0.16 0.1],'BackgroundColor',f_color,...
    'fontname','Times','FontSize',16,'Visible','off','ForegroundColor','black',...
    'HorizontalAlignment','center');

handles.EMGyes_button = uicontrol('Style','pushbutton','String','Yes',...
    'units','normalized','Position', [0.42 0.53 0.06 0.05],'Callback',@EMG_YES,...
    'fontname','Times','FontSize',16,'Visible','off');

handles.EMGno_button = uicontrol('Style','pushbutton','String','No',...
    'units','normalized','Position', [0.52 0.53 0.06 0.05],'Callback',@EMG_NO,...
    'fontname','Times','FontSize',16,'Visible','off');

% "Load movement" button
handles.loadMvmt_button = uicontrol('Style','pushbutton','String','Load Movement data',...
    'units','normalized','Position', [0.42 0.58 0.2 0.1],'Callback',@LoadMvmt,...
    'fontname','Times','FontSize',16,'Visible','off');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Block number and training state text

handles.blockText = uicontrol('Style','text','Visible','off','BackgroundColor',f_color,'String',...
    'Block:','ForegroundColor','black','FontSize',16,'units','normalized', 'Position',...
    [0.03 0.95 0.045 0.035],'Horizontalalignment','center','fontname','Times');

handles.BlockBox = uicontrol('Style','text','Visible','off','BackgroundColor',f_color,'String',...
    '','ForegroundColor','black','FontSize',21,'units','normalized', 'Position',...
    [0.03 0.905 0.045 0.035],'Horizontalalignment','center','fontname','Times');

handles.trainText = uicontrol('Style','text','Visible','off','BackgroundColor',f_color,'String',...
    'Training:','ForegroundColor','black','FontSize',16,'units','normalized', 'Position',...
    [0.09 0.95 0.055 0.035],'Horizontalalignment','center','fontname','Times');

handles.trainBox = uicontrol('Style','text','Visible','off','BackgroundColor',f_color,'String',...
    'Yes','ForegroundColor','black','FontSize',20,'units','normalized', 'Position',...
    [0.09 0.905 0.055 0.035],'Horizontalalignment','center','fontname','Times');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% block/resume training buttons

handles.stopTrain = uicontrol('Style','pushbutton','String','Stop training',...
    'units','normalized','Position', [0.16 0.945 0.10 0.045],'Callback',@stopTraining,...
    'fontname','Times','fontweight','bold','FontSize',14,'Visible','off',...
    'foregroundcolor',[237, 45, 50]./255);

handles.resumeTrain = uicontrol('Style','pushbutton','String','Resume training',...
    'units','normalized','Position', [0.16 0.895 0.10 0.045],'Callback',@resumeTraining,...
    'fontname','Times','fontweight','bold','FontSize',14,'Visible','off',...
    'foregroundcolor',[37, 191, 37]./255);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ML algorithm text and edit box
x_algo_text = .265;
w_algo_text = .08;

handles.mlAlgoText = uicontrol('Style','text','visible','off','Backgroundcolor',f_color,...
    'String','ML algorithm:','units','normalized','Position',[x_algo_text 0.948 w_algo_text 0.035],...
    'Horizontalalignment','right','FontSize',14,'fontname','Times');

handles.mlAlgoBox = uicontrol('Style','popupmenu','Visible','off','BackgroundColor',...
    'w','String',{'RF','ECOC'},'FontSize',14, 'units','normalized','Position',...
    [x_algo_text + w_algo_text, 0.946, w_algo_text, 0.04],...
    'Horizontalalignment','left','Callback',@changeAlgorithm);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% OOB error text and display box
x_oob_text = .265;
w_oob_text = .08;

handles.oobText = uicontrol('Style','text','visible','off','Backgroundcolor',f_color,...
    'String','OOB error:','units','normalized','Position',[x_oob_text 0.898 w_oob_text 0.035],...
    'Horizontalalignment','right','FontSize',14,'fontname','Times');

handles.oobBox = uicontrol('Style','text','Visible','off','BackgroundColor',...
    'w','String','','FontSize',16, 'units','normalized','Position',...
    [x_oob_text + w_oob_text + .006, 0.90, .068, 0.035],...
    'Horizontalalignment','center','fontname','Times');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% check video button

handles.checkvideo_button = uicontrol('Style','pushbutton','String','Check video',...
    'units','normalized','Position',[0.43 0.945 0.15 0.045],'Callback',@checkvidplayback,...
    'fontname','Times','Fontsize',16,'Visible','off');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% change state controls

handles.changeState = uicontrol('Style','pushbutton','String',sprintf('Change state'),...
    'units','normalized','Position', [0.43 0.895 0.15 0.045],'Callback',@changeState,...
    'fontname','Times','FontSize',17,'Visible','off');

handles.changeStateinput = uicontrol('Style','edit','Visible','off','BackgroundColor',...
    'w','String','','FontSize',14, 'units','normalized','Position',...
    [0.685 0.90 0.032 0.037],'Horizontalalignment','center','enable','off');

handles.statetext = uicontrol('Style','text','visible','off','Backgroundcolor',f_color,...
    'String','New state code:','units','normalized','Position',[0.58 0.885 0.1 0.05],...
    'Horizontalalignment','right','FontSize',15,'fontname','Times');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% go to block controls

handles.gotoblockText = uicontrol('Style','text','visible','off','Backgroundcolor',f_color,...
    'String','Go to block:','units','normalized','Position',[0.58 0.93 0.1 0.05],...
    'Horizontalalignment','right','FontSize',15,'fontname','Times');

handles.gotoblockBox = uicontrol('Style','edit','Visible','off','BackgroundColor',...
    'w','String','','FontSize',14, 'units','normalized','Position',...
    [0.685 0.945 0.032 0.037],'Horizontalalignment','center','Callback',@GoToBlock);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% next block controls

handles.NextBlock = uicontrol('Style','pushbutton','String','Next Block',...
    'units','normalized','Position', [0.725 0.895 0.16,0.095],'Callback',@nextblock,...
    'fontname','Times','FontSize',18,'Visible','off','fontweight','bold');


handles.NextBlock_noTrain = uicontrol('Style','pushbutton','String','<html>Next Block<br>(don''t train)',...
    'units','normalized','Position', [0.89 0.895 0.08 0.095],'Callback',@nextblock_notrain,...
    'fontname','Times','FontSize',15,'Visible','off');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% playback text

handles.playbackText = uicontrol('Style','text','Visible','on','BackgroundColor',f_color,'String',...
    '','ForegroundColor','black','FontSize',16,'units','normalized', 'Position',...
    [0.15 0.001 0.7 0.05],'Horizontalalignment','center','fontname','Times');

handles.GMTplaybackText = uicontrol('Style','text','Visible','on','BackgroundColor',f_color,'String',...
    '','ForegroundColor','black','FontSize',14,'units','normalized', 'Position',...
    [0.2 0.1 0.6 0.4],'Horizontalalignment','left','fontname','Times');

handles.GMTask_text = uicontrol('Style','text','String','Do GMT correction on Mvt data?',...
    'units','normalized','Position', [0.42 0.7 0.16 0.1],'BackgroundColor',f_color,...
    'fontname','Times','FontSize',16,'Visible','off','ForegroundColor','black',...
    'HorizontalAlignment','center');

handles.GMT_ask_no = uicontrol('Style','pushbutton','String','No','Visible','off',...
    'units','normalized','Position',[0.52 0.53 0.06 0.05],'Callback',@GMT_NO);

handles.GMT_ask_yes = uicontrol('Style','pushbutton','String','Yes','Visible','off',...
    'units','normalized','Position',[0.42 0.53 0.06 0.05],'Callback',@GMT_YES);



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% playback text

logo_axes = axes('position',[.9 .01 .085 .085]);
logo_image = image(mb_logo(:,:,1:3));
logo_image.AlphaData = 0.4;
set(logo_axes,'box','off','xtick',[],'xticklabel',[],...
    'ytick',[],'yticklabel',[],'Ycolor',f_color,'Xcolor',f_color,...
    'color',f_color);

%% CALLBACK FUNCTION DEFINITIONS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function LoadPrior(~,eventdata)
        % If the user has already done some coding, load the previous work and
        % make a note that this exists so as to (later) start at the correct
        % point in time.
        
        if strcmp(eventdata.Source.String,'Yes')
            % user input for the STATETIMES file
            [sFile, sDir] = uigetfile('.mat','Pick your previous statetimes.');
            
            % find number of file
            underscore_split = regexp(sFile,'_','split');
            handles.D.animal = underscore_split{1};
            matsplit = regexp(underscore_split{end},'\.','split');
            file_ID = matsplit{1};
            
            % find corresponding FEATURES file
            fFile = [handles.D.animal '_FEATURES_' file_ID '.mat'];
            
            % find the folder they are in
            statedat = load([sDir sFile]);
            featdat = load([sDir fFile]);
            
            % load the previously coded data
            handles.D.statetimes = statedat.statetimes;
            handles.D.training_data = featdat.training_set;
            fprintf('Loaded training data. Size: %u rows by %u cols.\n',...
                size(handles.D.training_data,1),size(handles.D.training_data,2));
            handles.D.prior = 1;
            handles.D.no_training = 0;
            
        else
            handles.D.prior = 0;
            handles.D.no_training = 0;
        end
        
        % set GUI for next step
        set(handles.loadLFP_button,'visible','on');
        set([handles.loadprioryes handles.loadpriorno handles.priortext],'visible','off');
        assignin('base','handles',handles)
        
    end


    function LoadLFPdata(~,~)
        
        % Ask user to select LFP file and load the data
        set(handles.loadLFP_button,'visible','off');
        [eFile, eDir] = uigetfile('','Select the LFP data file');
        set(handles.playbackText,'String','Loading LFP data. This will take a few minutes.','visible','on');
        drawnow
        temp = load([eDir eFile]);
        
        % set up EMG buttons
        set(handles.playbackText,'String','','visible','off');
        set(handles.useEMG_text,'visible','on');
        set(handles.EMGyes_button,'visible','on');
        set(handles.EMGno_button,'visible','on');
        
        % save LFP data in global handles variable
        handles.D.LFPdata = temp.LFPinfo;
        LFPstarts = [handles.D.LFPdata.startTime];
        
        % if no prior data, start at block 1
        if handles.D.prior == 0
            handles.D.block = 1; % CHANGE BACK TO 1
        else
            % otherwise find the correct block to start on
            a = [handles.D.LFPdata.startTime]';
            b = find(a>handles.D.statetimes(end,2),1);
            handles.D.block = b;
        end
        % set up other global parameters
        handles.D.AutoStates = [];
        handles.littlehorse = 0;
        handles.master_doTrain = 1;
        
        assignin('base','handles',handles);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
    end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % if no EMG
    function EMG_NO(~,~)
        
        set(handles.useEMG_text,'visible','off');
        set(handles.EMGyes_button,'visible','off');
        set(handles.EMGno_button,'visible','off');
        
        % set up NO_EMG flag
        handles.use_EMG = 0;
        
        set(handles.loadMvmt_button,'visible','on');
        assignin('base','handles',handles);
    end

    % if have EMG
    function EMG_YES(~,~)
        
        set(handles.useEMG_text,'visible','off');
        set(handles.EMGyes_button,'visible','off');
        set(handles.EMGno_button,'visible','off');
        
        % setup YES_EMG flag
        handles.use_EMG = 1;
        
        %  Index the EMG information by rec start times - - - - - - - - - -
        [handles.D.eDir]    = uigetdir('','Select folder with EMG data');
        % find animal name from this
        underscore_split  = regexp(handles.D.eDir,'_','split');
        animal_ID = underscore_split{end};
        % USER CHECK
        handles.D.animal = animal_ID;
        %handles.D.EMGcont   = dir(handles.D.eDir);
        handles.D.enames = dir([handles.D.eDir filesep '*EMG*']);
        handles.D.enames = {handles.D.enames.name};
        
        % find EMG data file names. Because the names have integer numbers
        % attached without leading 0s, they have to be indexed
        EMGindex = zeros(size(handles.D.enames,2),2);
        for ee = 1:size(handles.D.enames,2);
            
            t1 = []; t2 = [];
            
            EMGindex(ee,1) = ee;
            t1 = strfind(handles.D.enames{ee},'_EMG');
            t2 = strfind(handles.D.enames{ee},'.mat');
            EMGindex(ee,2) = str2double(handles.D.enames{ee}(t1+4:t2-1));
            
        end
        
        % sort the indexes
        EMGindex = sortrows(EMGindex,2);
        
        % store index to EMG files here
        handles.D.EMGidx = EMGindex;
        
        % set up load Mvt buttons
        set(handles.loadMvmt_button,'visible','on');
        assignin('base','handles',handles);
    end

    % if NO GMT correction
    function GMT_NO(~,~)
        
        set(handles.GMTask_text,'visible','off');
        set(handles.GMT_ask_yes,'visible','off');
        set(handles.GMT_ask_no,'visible','off');
        
        handles.do_GMT_correct = 0;
        
        assignin('base','handles',handles);
        uiresume(gcbf);
    end
    
    % if YES GMT correction
    function GMT_YES(~,~)
        
        set(handles.GMTask_text,'visible','off');
        set(handles.GMT_ask_yes,'visible','off');
        set(handles.GMT_ask_no,'visible','off');
        
        handles.do_GMT_correct = 1;
        
        assignin('base','handles',handles);
        uiresume(gcbf);
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Load movement
    function LoadMvmt(~,~)
        
        %  Load automated movement analysis - - - - - - - - - - - - - - - -
        set(handles.loadMvmt_button,'visible','off');
        % Ask user to find movement data file
        [eFile, eDir] = uigetfile('','Select movement data');
        
        % load the movement data file
        temp = load([eDir eFile]);
        if isfield(temp,'DATA')
            
            % load the smoothed movement and head movement data
            pmove = temp.DATA.smooth_movement';
            quiet = temp.DATA.quiet_mvt';
            
            % clean the data
            try
                handles.D.movement  = [pmove temp.DATA.frame_times(1:end)];
                handles.D.quiet  = [quiet temp.DATA.frame_times(1:end)];
                
%                 keyboard;
                badz = find(handles.D.movement(:,2) < handles.D.movement(1,2));
                handles.D.movement(badz,:) = [];
                handles.D.quiet(badz,:) = [];
            catch
                keyboard
            end
            
            if isfield(temp.DATA,'mask')
                handles.D.ROI       = temp.DATA.mask;
            else
                handles.D.ROI = [];
            end
            
            % check for GMT time.
            lfp_start = handles.D.LFPdata(1).startTime;
            lfp_end = handles.D.LFPdata(end).startTime;
            mvt_start = handles.D.movement(1,2);
            mvt_end = handles.D.movement(end,2);
            
            % display warning.
            % This happens because the LFP data is usually in local time,
            % while the movement data comes out in GMT time. They are
            % usually 5 hours apart, although for 2 weeks out of the year
            % this difference can be 4 (or 6?) hours. The code below fixes
            % the GMT offset no matter what time of year it is.
            warn_Str = sprintf(['WARNING! It looks like the LFP times and movement data',...
                ' times are misaligned. Possibly due to missing GMT time correction.\n',...
                'The start times are:\n%s for LFP\n%s for movement.\n\n',...
                'The end times are:\n%s for LFP\n%s for movement.\n\n',...
                '\n_____     *** Fixing GMT offset ***     _____\n\n\n'],...
                datestr(unixtime(lfp_start)),datestr(unixtime(mvt_start)),...
                datestr(unixtime(lfp_end)),datestr(unixtime(mvt_end)));
            set(handles.GMTplaybackText,'visible','on','string',warn_Str);
            set(handles.GMT_ask_yes,'visible','on');
            set(handles.GMT_ask_no,'visible','on');
            set(handles.GMTask_text,'visible','on');
            
            uiwait(gcf);
            
            if handles.do_GMT_correct
                
                % Movement timestamps are in GMT time. Adjust this to match
                % LFP timestamps.
                first_timestamp = unixtime(handles.D.movement(1,2));
                in_local_time = TimezoneConvert_ATP(first_timestamp,'UTC','America/New_York');
                in_local_unixtime = unixtime(datevec(in_local_time));
                hour_offset = round(abs(handles.D.movement(1,2) - in_local_unixtime))/3600;
                handles.D.hour_offset = hour_offset;
                handles.D.movement(:,2) = handles.D.movement(:,2) - hour_offset*3600;
                handles.D.quiet(:,2) = handles.D.quiet(:,2) - hour_offset*3600;
                
                correct_str = sprintf(['The movement timestamps are off by %u hours.\n',...
                    'These are the new movement start and end times:\n',...
                    '%s (start);\n  %s (end)\n\n',...
                    'These are the LFP start and end times:\n',...
                    '%s (start);\n  %s (end)\n\n'],hour_offset,...
                    datestr(unixtime(handles.D.movement(1,2))),...
                    datestr(unixtime(handles.D.movement(end,2))),...
                    datestr(unixtime(lfp_start)),...
                    datestr(unixtime(lfp_end)));
                
                set(handles.GMTplaybackText,'string',correct_str);
                pause(2.5);
                set(handles.GMTplaybackText,'visible','off');
                
            end
        else
            handles.D.movement = [temp.outdata.DATA.smooth_movement temp.outdata.DATA.frame_times(2:end)];
            handles.D.quiet = [temp.outdata.DATA.quiet temp.outdata.DATA.frame_times(2:end)];
            try
                handles.D.ROI = temp.outdata.ROI;
            catch
                handles.D.ROI = temp.outdata.DATA.mask;
            end
        end
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        assignin('base','handles',handles);
        processData;
        

    end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%                         Finished UI Load Data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % process the data
    function processData(~,~)
        
        set(handles.playbackText,'string','Processing data.');
        drawnow;
        try
            cla(handles.LFPfig);
        end
        
        % initialize variables
        handles.D.semgY         = [];
        handles.D.emgData       = [];
        handles.D.P             = [];
        handles.D.F             = [];
        handles.D.pX            = [];
        handles.D.rs_emg        = [];
        handles.D.statetimes    = [];
        handles.D.firstdata     = 10;
        
        % check to make sure there's LFP data in this block. if not, find
        % the next occurence of LFP data.
        if isempty(handles.D.LFPdata(handles.D.block).spectrogram)
            dumb = 0;
            while dumb == 0
                handles.D.block = handles.D.block+1;
                
                if ~isempty(handles.D.LFPdata(handles.D.block).spectrogram)
                    dumb = 1;
                    handles.D.firstdata = handles.D.block;
                end
            end
        end
        
        set(handles.gotoblockBox,'string','?');
        
        set(handles.BlockBox,'visible','on','string',num2str(handles.D.block));
        set([handles.blockText handles.changeState handles.gotoblockText handles.gotoblockBox ...
            handles.checkvideo_button handles.mlAlgoBox handles.mlAlgoText ...
            handles.oobBox handles.oobText handles.statetext handles.changeStateinput],'visible','on');
        handles.D.ml_algorithm = handles.mlAlgoBox.String{handles.mlAlgoBox.Value};
        handles.mlAlgoBox.Enable = 'off';
        handles.changeStateinput.Enable = 'off';
        
        %% ATP NOTE 02
        % Load and process LFP to get power for delta and theta frequency
        % bands. NOTE: may use z-score here too, but these will have to be
        % calculated block-by-block. Might work, might not.
        
        % Load LFP data (needs to be loaded first)
        handles.D.P = handles.D.LFPdata(handles.D.block).spectrogram.P;
        handles.D.F = handles.D.LFPdata(handles.D.block).spectrogram.F;
        % these are the timestamps for the LFP data
        handles.D.pX = linspace(handles.D.LFPdata(handles.D.block).startTime,  handles.D.LFPdata(handles.D.block).startTime...
            + handles.D.LFPdata(handles.D.block).duration, size(handles.D.P,2) );
        
        %Old method:
        %{
        deltapower  = sum(10*log10(abs(handles.D.P(handles.D.F<=4,:))));
        deltapower  = deltapower/abs(max(deltapower));
        deltapower  = 10*(deltapower+abs(min(deltapower))); % make positive and scale for the plot
        deltapower  = deltapower.^1.5;
        deltapower  = smooth(deltapower,100);
        
        thetapower  = sum(10*log10(abs(handles.D.P(handles.D.F>=4 & handles.D.F<=8,:))));
        thetapower  = thetapower/abs(max(thetapower));
        thetapower  = 10*(thetapower+abs(min(thetapower)));
        thetapower  = smooth(thetapower,200);
        thetapower  = thetapower.^2;
        %}
        
        % Freq ranges:
        % DELTA: 0 - 4 Hz
        delta_freq = [0 4];
        % THETA: 5 - 8 Hz
        theta_freq = [5 8];
        % GAMMA: 40-100 Hz
        gamma_freq = [40 100];
        handles.use_gamma = 1;
        
        % New method:
        % Get absolute power from Power Spectral Density, for delta,
        % theta, and whole range.
        delta_absolute_power = sum(handles.D.P(handles.D.F >= delta_freq(1) & handles.D.F <= delta_freq(2),:));
        theta_absolute_power = sum(handles.D.P(handles.D.F >= theta_freq(1) & handles.D.F <= theta_freq(2),:));
        total_absolute_power = sum(handles.D.P(:,:));
        if handles.use_gamma
            gamma_absolute_power = sum(handles.D.P(handles.D.F >= gamma_freq(1) & handles.D.F <= gamma_freq(2),:));
        end
        
        % Get relative power fraction for each frequency band and smooth
        % delta
        delta_relative_power = delta_absolute_power ./ total_absolute_power;
        deltapower = smooth(delta_relative_power,200);
        delta_mean = nanmean(deltapower);
        delta_std  = nanstd(deltapower);
        delta_zscore = (deltapower - delta_mean) ./ delta_std;
        % theta
        theta_relative_power = theta_absolute_power ./ total_absolute_power;
        thetapower = smooth(theta_relative_power,200);
        theta_mean = nanmean(thetapower);
        theta_std  = nanstd(thetapower);
        theta_zscore = (thetapower - theta_mean) ./ theta_std;
        % gamma
        if handles.use_gamma
            gamma_relative_power = gamma_absolute_power ./ total_absolute_power;
            gammapower = smooth(gamma_relative_power,200);
            gamma_mean = nanmean(gammapower);
            gamma_std  = nanstd(gammapower);
            gamma_zscore = (gammapower - gamma_mean) ./ gamma_std;
        end
        
        % save relevant variables to handles
        handles.D.deltapower = deltapower;
        handles.D.thetapower = thetapower;
        handles.D.delta_z    = delta_zscore;
        handles.D.theta_z    = theta_zscore;
        if handles.use_gamma
            handles.D.gammapower = gammapower;
            handles.D.gamma_Z = gamma_zscore;
        end
        
        
        %% ATP NOTE 01
        % Load and process EMG signal.
        
        if handles.use_EMG
            % Load EMG data
            etemp = load([handles.D.eDir  filesep handles.D.enames{ handles.D.EMGidx(handles.D.block,1) } ]);
            emgLtemp = etemp.EMGdata.f_EMG_L;
            emgRtemp = etemp.EMGdata.f_EMG_R;
            if isempty(emgLtemp) || isempty(emgRtemp)
                have_both_EMGs = 0;
            else
                have_both_EMGs = 1;
            end
            fprintf('Have both EMGs: %u.\n',have_both_EMGs);
            % old method:
            %{
        if var(emgLtemp) > var(emgRtemp);
            handles.D.emgData = emgLtemp;
        else
            handles.D.emgData = emgRtemp;
        end
            %}
            
            use_EMG_avg = 0;
            
            
            if ~use_EMG_avg && have_both_EMGs
                % alternative to averaging EMG together: subtract one from
                % other
                emgavg = (emgLtemp - emgRtemp);
            else
                % Average together the signal from both hemispheres.
                emgavg = (emgRtemp + emgLtemp) / 2;
            end
            handles.D.emgData = emgavg;
            
            %         keyboard;
            
            % Smooth and resample EMG signal
            handles.D.emgX = linspace(handles.D.LFPdata(handles.D.block).startTime,  handles.D.LFPdata(handles.D.block).startTime...
                + handles.D.LFPdata(handles.D.block).duration, size(handles.D.emgData,2) );
            handles.D.emgY = handles.D.emgData;
            
            %         smooth_span = 200;
            %         handles.D.semgY = smooth(handles.D.emgY,smooth_span,'rloess');
            
            % take a moving average to resample the EMG data (to match LFP length)
            %         bb = round(linspace(1, size(handles.D.semgY,1), size(handles.D.pX,2)+1  ));
            
            bb = round(linspace(1, size(handles.D.emgY,2), size(handles.D.pX,2)+1  ));
            
            for ee = 1:size(bb,2)-1
                %             handles.D.rs_emg(ee) = nanmean( handles.D.semgY( bb(ee):bb(ee+1) , 1 ));
                handles.D.rs_emg(ee) = nanmean(handles.D.emgY( 1, bb(ee):bb(ee+1) ) );
            end
            
            rs_emg = handles.D.rs_emg;
            % Transform EMG to z-score for generalizability. NOTE: this is
            % z-score for this block. May need to extend to whole experiment to
            % make truly generalizable.
            emg_mean = nanmean(rs_emg);
            emg_std  = nanstd(rs_emg);
            emg_zscore = (rs_emg - emg_mean) ./ emg_std;
            handles.D.emg_z = emg_zscore;
            
        else
            fprintf('\t*** NOT USING EMG DATA! ***\n');
            
            emg_zscore = zeros(size(deltapower,1),1);
            rs_emg = zeros(size(deltapower,1),1);
            
            
        end
        
        
        
        %% ATP NOTE 03
        % The following used to scale theta and delta power to set their means
        % equal to 2. Since I am now using relative power fractions,
        % shouldn't need this anymore. Commented out.
        
        % Old method:
        %{
        % examine data from the first five hours of the recording to figure
        % out how to scale the LFP and EMG data from individual animals
        % such that the data will be approachable via a generalized code
        if isempty(handles.D.deltafactor);
            dfac = []; tfac = []; eefac = [];
            
            for ee = handles.D.firstdata:handles.D.firstdata+4;
                
                dtemp = []; P = []; F = []; emgtemp = [];
                
                P = handles.D.LFPdata(ee).spectrogram.P;
                F = handles.D.LFPdata(ee).spectrogram.F;
                
                dtemp  = sum(10*log10(abs(P(F<=4,:))));
                dtemp  = dtemp/abs(max(dtemp));
                dtemp  = 10*(dtemp+abs(min(dtemp))); % make positive and scale for the plot
                dtemp  = dtemp.^1.5;
                dtemp  = smooth(dtemp,100);
                
                ttemp  = sum(10*log10(abs(P(F>=6 & F<=8,:))));
                ttemp  = ttemp/abs(max(ttemp));
                ttemp  = 10*(ttemp+abs(min(ttemp)));
                ttemp  = smooth(ttemp,200);
                ttemp  = ttemp.^2;
                
                dfac    = [dfac; dtemp];
                tfac    = [tfac; ttemp];
                
            end
            
            handles.D.deltafactor = mean(dfac)/2;   % Set the mean delta power to 2
            handles.D.thetafactor = mean(tfac)/2;
        end
        
        deltapower  = deltapower/handles.D.deltafactor;
        thetapower  = thetapower/handles.D.thetafactor;
        %}
        
        %% ATP NOTE 04
        % Previously we scaled the EMG for display purposes. Old code here.
        
        % Old method:
        %{
        emgX = linspace(handles.D.LFPdata(handles.D.block).startTime,  handles.D.LFPdata(handles.D.block).startTime...
            + handles.D.LFPdata(handles.D.block).duration, size(handles.D.emgData,2) );
        emgY = 10.*(handles.D.emgData/max(handles.D.emgData));
        
        handles.D.semgY = smooth(handles.D.emgData,2000);
        handles.D.semgY = 10.*(handles.D.semgY/max(handles.D.semgY));
        handles.D.semgY = handles.D.semgY/mean(handles.D.semgY);
        
        % take a moving average to smooth/resample the EMG
        bb = round(linspace(1, size(handles.D.semgY,1), size(handles.D.pX,2)+1  ));
        
        for ee = 1:size(bb,2)-1;
            handles.D.rs_emg(ee) = mean( handles.D.semgY( bb(ee):bb(ee+1) ));
        end
        %}
        
        
        
        %% ATP NOTE 05
        % Processing movement data. For more generalizability, use z-score
        % of movement trace.
        
        % Old method:
        %{
        % Process the video analysis data for plotting:
        mvmt = handles.D.movement;
        mvmt(:,1) = mvmt(:,1) ./ (0.5*nanmean( mvmt( mvmt(:,1)~=0 ,1)   ) ); % DOUBLING OUTPUT WITH 0.5 TO ADJUST FOR ATP NEW CODE 7/27/16 KBH
        mvmt(mvmt(:,2)<handles.D.pX(1),:)   = [];
        mvmt(mvmt(:,2)>handles.D.pX(end),:) = [];
%         mvmt(:,1) = mvmt(:,1)/mean( mvmt( mvmt(:,1)~=0 ,1)   );
        %}
        
        % clear zscore
        mvt_zscore = [];
        
        % get movement trace for whole experiment
        mvmt = handles.D.movement;
        
        
        mvt_mean = nanmean(mvmt(:,1));
        mvt_std  = nanstd(mvmt(:,1));
        
        mvt_zscore = (mvmt(:,1) - mvt_mean) ./ mvt_std;
        
        % remove times outside of current block
%         keyboard;
        prev_times = find(mvmt(:,2) < handles.D.pX(1));
        post_times = find(mvmt(:,2) > handles.D.pX(end));
        mvt_zscore([prev_times; post_times],:)   = [];
        %         mvt_zscore(mvmt(:,2)>handles.D.pX(end),:) = [];
        
        %         bb = round(linspace(1, size(mvt_zscore,1), size(handles.D.pX,2)+1  ));
        
        % keyboard;
        
        % check for empty movement arrays. this can happen if the video cut
        % out for example
        
        if isempty(prev_times)
            rs_mvt_z = nan(size(handles.D.deltapower));
        else
            if isempty(post_times)
                post_times = length(mvmt);
            end
            mvt_diff = diff(mvmt(prev_times(end)+1:post_times(1)-1,2));
            if numel(find(mvt_diff > 60)) > 2
                disp('Found multiple skipped timestamps in the movement!!!');
            end
            
            [~,maxidx] = max(mvt_diff);
            if isempty(maxidx)
                rs_mvt_z = nan(1, length(mvt_diff));
            else
                px_idx = find(handles.D.pX(1,:) > mvmt(prev_times(end)+maxidx,2),1,'first');
                px_idx2 = find(handles.D.pX(1,:) < mvmt(prev_times(end)+maxidx+1,2),1,'last');
                bb = zeros(1,size(handles.D.pX,2)+1);
                bb(1:px_idx) = linspace(1,maxidx,px_idx);
                bb(px_idx+1:px_idx2) = bb(px_idx);
                bb(px_idx2+1:end) = linspace(maxidx+1,size(mvt_zscore,1),size(bb,2)-px_idx2);
                
                for ee = 1:size(bb,2)-1
                    try %THIS WAS SOPHIE
                        rs_mvt_z(ee) = nanmean( mvt_zscore( bb(ee):bb(ee+1) , 1));
                    catch
                        keyboard
                    end
                end
                
                rs_mvt_z(px_idx+1:px_idx2) = nan(1,px_idx2-px_idx);
            end
            
            % adjust so that no movement = 0
            rs_mvt_z = rs_mvt_z + abs(min(rs_mvt_z));
        end
        handles.D.rs_mvt_z = rs_mvt_z;
        
        % Do same thing for quiet (head only) movement ~ Sam Wacks
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        % clear zscore
        quiet_zscore = [];
        % get movement trace for whole experiment
        quiet = handles.D.quiet;
        quiet_mean = nanmean(quiet(:,1));
        quiet_std  = nanstd(quiet(:,1));
        quiet_zscore = (quiet(:,1) - quiet_mean) ./ quiet_std;
        % remove unwanted times
        prev_times = find(quiet(:,2) < handles.D.pX(1));
        post_times = find(quiet(:,2) > handles.D.pX(end));
        quiet_zscore([prev_times; post_times],:)   = [];
        
        if isempty(prev_times)
            rs_quiet_z = nan(size(handles.D.deltapower));
        else
            if isempty(post_times)
                post_times = length(quiet);
            end
            quiet_diff = diff(quiet(prev_times(end)+1:post_times(1)-1,2));
            if numel(find(quiet_diff > 60)) > 2
                disp('Found multiple skipped timestamps in the movement!!!');
            end
            
            [~,maxidx] = max(quiet_diff);
            if isempty(maxidx)
                rs_quiet_z = nan(1, length(quiet_diff));
            else
                px_idx = find(handles.D.pX(1,:) > quiet(prev_times(end)+maxidx,2),1,'first');
                px_idx2 = find(handles.D.pX(1,:) < quiet(prev_times(end)+maxidx+1,2),1,'last');
                
                bb = zeros(1,size(handles.D.pX,2)+1);
                bb(1:px_idx) = linspace(1,maxidx,px_idx);
                bb(px_idx+1:px_idx2) = bb(px_idx);
                bb(px_idx2+1:end) = linspace(maxidx+1,size(quiet_zscore,1),size(bb,2)-px_idx2);
                
                for ee = 1:size(bb,2)-1
                    try %THIS WAS SOPHIE
                        rs_quiet_z(ee) = nanmean(quiet_zscore( bb(ee):bb(ee+1) , 1));
                    catch
                        keyboard
                    end
                end
                
                rs_quiet_z(px_idx+1:px_idx2) = nan(1,px_idx2-px_idx);
            end
        end
        % adjust so that no movement = 0
        %rs_quiet_z = rs_quiet_z + abs(min(rs_quiet_z));
        handles.D.rs_quiet_z = rs_quiet_z;
        
        
        %% ATP NOTE 06
        % This section calculates the features for the classifier model.
        %
        % FEATURES TO BE USED
        % _ EMG
        % _ Delta
        % _ Theta
        % _ Mvt
        % _ (Delta - Theta)
        % _ var(EMG)
        % _ var(Mvt)
        % _ quiet mvt
        
        % Theta-Delta difference
        dt_diff = deltapower - thetapower;
        handles.D.dt_diff = dt_diff;
        
        % make sure arrays are in similar orientation
        if size(rs_mvt_z,1) ~= size(deltapower,1) && size(rs_mvt_z,1) == size(deltapower,2)
            rs_mvt_z = rs_mvt_z';
        end
        if size(rs_quiet_z,1) ~= size(deltapower,1) && size(rs_quiet_z,1) == size(deltapower,2)
            rs_quiet_z = rs_quiet_z';
        end
        if size(thetapower,1) ~= size(deltapower,1) && size(thetapower,1) == size(deltapower,2)
            thetapower = thetapower';
        end
        if size(emg_zscore,1) ~= size(deltapower,1) && size(emg_zscore,1) == size(deltapower,2)
            emg_zscore = emg_zscore';
        end
        if size(rs_emg,1) ~= size(deltapower,1) && size(rs_emg,1) == size(deltapower,2)
            rs_emg = rs_emg';
        end
        if size(dt_diff,1) ~= size(deltapower,1) && size(dt_diff,1) == size(deltapower,2)
            dt_diff = dt_diff';
        end
        
        % COMPILE DATA REQUIRED FOR FEATURE CALCULATION IN A MATRIX
        handles.D.feature_data = [deltapower, thetapower, dt_diff, rs_mvt_z, rs_quiet_z];
        if handles.use_EMG
            handles.D.feature_data = [handles.D.feature_data, rs_emg];
        end
        if handles.use_gamma
            handles.D.feature_data = [handles.D.feature_data, gammapower, rs_emg];
        end
        
        
        % Old method:
        %{
        % put the relevant data into a single matrix for easy passing between
        % functions
        sm_delta = smooth(deltapower,20);
        handles.D.seq = [sm_delta thetapower handles.D.rs_emg'];
        %}
        
        %% ATP NOTE 08
        % Depending on block number, either run tensectest or ML coding
        
        % On first 10 (?) blocks, just use tensectest_new (new version)
        handles.D.X_Block = 10;
        
        % if we are not yet at the RF algorithm section, use a previous
        % threshold-based model to estimate states. This is basically
        % useless, ends up being almost completely manual coding
        testdat = [];
        if handles.D.block <= handles.D.X_Block
            sampling_rate = 2; % Hz
            set(handles.playbackText,'visible','on','string','Scoring data');
            [testdat, rawdat, block_features] = score_10sec_noMLS_Quiet(handles.D.feature_data, handles.D.pX, sampling_rate);
            handles.D.features{handles.D.block} = block_features;
        else
            if handles.D.block > handles.D.X_Block
                set([handles.stopTrain, handles.resumeTrain, ...
                    handles.trainBox, handles.trainText handles.playbackText],'visible','on');
                handles.mlAlgoBox.Enable = 'on';
            end
            % ALGORITHM HERE
            % handles.D.ml_algorithm specifies which algorithm to use.
            %  - 'RF' : Random Forest
            %  - 'ECOC' : Error-correcting output code - THIS IS NOT
            %             IMPLEMENTED YET
            
            if handles.master_doTrain
                handles.D.trained_Mdl_path = [];
                
                % train model based on previous blocks
                
                set(handles.playbackText,'string',sprintf('Training and saving %s model...',handles.D.ml_algorithm));
                drawnow;
                talg0 = tic;
                RF_nTrees = [];
                [pymdl_path,oob_err,savemdl_dir] = AutoVidCode_trainModel_PyWrap(handles.D.training_data,...
                    handles.D.animal,handles.D.ml_algorithm,RF_nTrees,myname);
                set(handles.oobBox,'String',sprintf('%.2f%%',100*(1-oob_err)));
                handles.OOB_ERROR(handles.D.block) = (1-oob_err)*100;
                
                %--------------------------------------------------------------
                % ******________________ IMPORTANT NOTE: ________________******
                % By default the above function uses n=200 trees. To change
                % the number of trees, pass it as the 4th argument to the
                % function.
                %--------------------------------------------------------------
                %--------------------------------------------------------------
                
                % THIS SHOULD ALWAYS BE 1!!! 0 ONLY FOR DEBUGGING
                do_save_mdl = 1;
                if do_save_mdl
                    Mdl_info.ML_algo = handles.D.ml_algorithm;
                    Mdl_info.block = handles.D.block;
                    
                    % Saving model
                    fprintf('\n  Saving your ML model...\n');
                    ts0 = tic;
                    if ~exist(savemdl_dir), mkdir(savemdl_dir); end
                    save([savemdl_dir filesep 'pyMdl_info.mat'],'Mdl_info');
                    copyfile(pymdl_path,savemdl_dir);
                    ts1 = toc(ts0);
                    fprintf('Done! Time elapsed: %.2f seconds.\n\n',ts1);
                end
                
                talg1 = toc(talg0);
                set(handles.playbackText,'string',sprintf('Training and saving %s model... That took %.2f seconds.',...
                    handles.D.ml_algorithm,talg1));
                handles.D.trained_Mdl_path = pymdl_path;
            else
                set(handles.playbackText,'string',sprintf('Training is OFF! Using last trained model.'));
            end
            pause(.1);
            
            %             loaded_Mdl = AutoVidCode_loadModel(handles.D.animal,ml_algorithm);
            
            % parameters for video scoring. Bin size for scoring and
            % sampling rate of input data
            bin_sz = 10; % seconds
            sampling_rate = 1; % Hz
            set(handles.playbackText,'visible','on','string',...
                sprintf('Auto-scoring using %s model. Please wait...',handles.D.ml_algorithm));
            % Call to autoscoring function
            [testdat, rawdat, block_features] = auto_video_score_PyWrap_Quiet(handles.D.feature_data, handles.D.trained_Mdl_path, ...
                handles.D.pX, bin_sz, sampling_rate, myname);
            handles.D.features{handles.D.block} = block_features;
        end
        
        
        
        % Old method:
        %{
        % A more generalizable way to do the 10-sec test would be to look
        % at variations from the mean (in each block or across blocks) to
        % determine thresholds for active/quiet waking in MVT or for NREM
        % in delta or REM in theta, etc.
        
        % TRY RUNNING WITH ONLY THIS LINE AS VERSION 2:
        testdat = [];
        testdat = tensectest(handles.D.seq, handles.D.pX, mvmt);
        
        tD      = diff(testdat(:,1));
        tD(end) = 1;
        killem  = find(tD == 0);
        killem  = killem+1;
        testdat2 = testdat;
        testdat2(killem,:) = [];
        %}
        
        
        
        %%
        % - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - -
        %  - - - - - - - - - - - DO THE PLOTTING - - - - - - - - - - - - -
        % - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - -
        
        % CLEAR AXES
        %         keyboard;
        try
            tmph = get(handles.f,'children');
            for uu=1:max(size(tmph))
                if strcmp(tmph(uu).Type,'axes')
                    cla(tmph(uu));
                    delete(tmph(uu));
                end
            end
            if exist('handles.LFPfig','var')
                delete(handles.LFPfig);
                if handles.use_EMG
                    delete(handles.EMGfig);
                end
                delete(handles.feat_axes);
            end
            drawnow;
        catch
            keyboard;
        end
        
        % Plot the spectrogram:
        handles.LFPfig = axes('Position',[.05 .13 .9 .75],'visible','off','Parent',handles.f);
        colormap parula
        
        axes(handles.LFPfig);
        psd_axes_pos = handles.LFPfig.Position;
        max_LFP_y = 148; % 148 = 15 Hz
        imagesc(handles.D.pX,( handles.D.F(1:max_LFP_y) ),10*log10(abs(handles.D.P(1:max_LFP_y,:))));
        set(handles.LFPfig,'ydir','normal','visible','off')
        %xtemp = (handles.D.pX - (handles.D.pX(1)))/60;
        
        
        
        % old:
        %{
        plot(handles.D.pX, thetapower,'k','linewidth',2);
        plot(handles.D.pX, sm_delta,'r','linewidth',2);
        plot(mvmt(:,2),mvmt(:,1),'m--','linewidth',2);
        plot(handles.D.pX, handles.D.rs_emg,'.','linewidth',0.5);
        %}
        
        % plot the z-score of movement and EMG on right y axis
        %         keyboard;
        if handles.use_EMG
            handles.EMGfig = axes('Position',[.05 .13 .9 .3],'visible','off','Parent',handles.f);
            axes(handles.EMGfig);
            emg_ylim = abs(mean(handles.D.rs_emg)) + 0.5*std(handles.D.rs_emg);
            yscale_emg = 'log';
            if min(handles.D.rs_emg) < 0
                emg_ylim_1 = mean(handles.D.rs_emg) - 3*std(handles.D.rs_emg);
                emg_ylim_2 = mean(handles.D.rs_emg) + 3*std(handles.D.rs_emg);
                yscale_emg = 'linear';
            else
                emg_ylim_1 = -emg_ylim;
                emg_ylim_2 = emg_ylim;
            end
            plot_emg = plot(handles.D.rs_emg,'-','linewidth',2,'color',[0 .45 .74]);
            set(handles.EMGfig,'xlim',[0 size(handles.D.rs_emg,2)],'xtick',[],'xticklabel',[],...
                'yscale',yscale_emg,...);
                'ylim',[emg_ylim_1 emg_ylim_2]);
            set(handles.EMGfig,'visible','off');
        end
        
        handles.feat_axes = axes('Position',psd_axes_pos,'Color','none');
        yyaxis right; hold on;
        pm = false;
        if length(handles.D.pX) == length(handles.D.rs_quiet_z)
            pm = true;
            plot_quiet = plot(handles.D.pX, handles.D.rs_quiet_z, 'g--', 'linewidth', 1);
        end
        if length(handles.D.pX) == length(handles.D.rs_mvt_z)
            plot_mvt = plot(handles.D.pX, handles.D.rs_mvt_z,'m--','linewidth',2);
        end
        
        %         plot_emg = plot(handles.D.pX, handles.D.emg_z,'.','linewidth',0.5,...
        %             'color',[0 .45 .74]);
        hold off;
        % plot fraction of delta and theta power on left y axis (and
        % delta-theta)
        yyaxis left; hold on;
        set(gca,'visible','on','fontsize',14);
        plot_theta = plot(handles.D.pX, handles.D.thetapower,'-k','linewidth',2);
        plot_delta = plot(handles.D.pX, handles.D.deltapower,'-r','linewidth',2);
        if handles.use_gamma
            plot_gamma = plot(handles.D.pX,handles.D.gammapower,'-y','linewidth',2);
        end
        plot_dt = plot(handles.D.pX, handles.D.dt_diff,'--c','linewidth',1.5);
        hold off;
        
        % legend
        if pm
            if handles.use_EMG
                if handles.use_gamma
                    plot_legend = legend([plot_theta,plot_delta,plot_gamma,plot_mvt,plot_quiet,plot_emg],...
                        {'Theta','Delta','Gamma','Movement','Quiet Movement','EMG'},'fontsize',14);
                else
                    plot_legend = legend([plot_theta,plot_delta,plot_mvt,plot_quiet,plot_emg],...
                        {'Theta','Delta','Movement','Quiet Movement','EMG'},'fontsize',14);
                end
            else
                plot_legend = legend([plot_theta,plot_delta,plot_mvt,plot_quiet],...
                    {'Theta','Delta','Movement','Quiet Movement'},'fontsize',14);
            end
        else
            if handles.use_EMG
                if handles.use_gamma
                    plot_legend = legend([plot_theta,plot_delta,plot_gamma,plot_emg],...
                        {'Theta','Delta','Gamma','EMG'},'fontsize',14);
                else
                    plot_legend = legend([plot_theta,plot_delta,plot_emg],...
                        {'Theta','Delta','EMG'},'fontsize',14);
                end
            else
                plot_legend = legend([plot_theta,plot_delta,],...
                    {'Theta','Delta'},'fontsize',14);
            end
        end
        set(plot_legend,'position',[0.85 0.66 0.08 0.07],'color',[1 1 1]);
        
        % set the x axis labels to minutes elapsed in the plotted dataframe
        xlims        = get(handles.feat_axes,'xlim');
        xlims2 = get(handles.LFPfig,'xlim');
        onesec_inpts = size(handles.D.pX,2)/3600;
        set(handles.feat_axes,'xlim',xlims2);
        set(handles.LFPfig,'xlim',xlims2);
        if handles.use_EMG
            set(handles.EMGfig,'xtick',[]);
        end
        set(handles.feat_axes,'XTick',[xlims2(1):600:xlims2(2)]);
        currlabel   = get(handles.feat_axes,'XTick');
        newlabel    = currlabel - xlims2(1);
        set(gca,'XTickLabel',newlabel/60);
        
        xl          = xlabel('Time (minutes)');
        yyaxis left
        yl_1        = ylabel('Fraction of PSD (LFP features)');
        yyaxis right
        yl_2        = ylabel('Z-score (Movement feature)');
        set([xl yl_1 yl_2],'Fontname','Myriad','Fontsize',16)
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Plot the algorithmically coded video:
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Ay = 8;
        my_ymin = -2;
        axes(handles.feat_axes);
        
        % colors
        clr = [0.5 0.5 0.5; 0.2 0.2 0.3; NaN NaN NaN; 0.4 1.0 0.2; 1 1 0];
        
        % plot the estimated state codes
        for ee = 1:size(testdat,1)
            t0      = testdat(ee,2);
            t1      = testdat(ee,3);
            
            try
                yyaxis right
                rectangle( 'Position', [t0, Ay, t1 - t0 , 1],...
                    'facecolor',clr(testdat(ee),:),'linestyle','none');
            catch
                disp('Caught in plotting auto output');
                keyboard
            end
        end
        
        handles.D.rawdat = testdat;
        
        yyaxis right
        set(gca,'ylim',[my_ymin Ay+1]);
        yyaxis left
        set(gca,'ylim',[0 1]);
        set(gca,'TickLength',[0 0.025]);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        % Eliminate repeated bins:
        stdiff = [1; diff(testdat(:,1))];
        testdat(stdiff == 0,:) = [];
        
        % Check to see if the user is recoding this block of data - if so,
        % clear the old entries and write in the new data:
        sortflag = 0;
        if ~isempty(handles.D.AutoStates)
            
            binS = handles.D.pX(1);
            binE = handles.D.pX(end);
            
            if any(handles.D.AutoStates(:,2) >=binS & handles.D.AutoStates(:,2) <binE)
                handles.D.AutoStates(handles.D.AutoStates(:,2) >=binS & handles.D.AutoStates(:,2) <binE,:) = [];
                sortflag = 1;
            end
        end
        
        % Now write output to handles
        handles.D.AutoStates = [handles.D.AutoStates; testdat(:,1:2)];
        
        if sortflag == 1
            handles.D.AutoStates = sortrows(handles.D.AutoStates,2);
        end
        
        % figure out whether to save or move on
        
        set([handles.NextBlock, handles.NextBlock_noTrain],'visible','on');
        
        set(handles.playbackText,'string',[]);
        
        assignin('base','handles',handles);
        
        %% ATP NOTE 11
        % Before end of function, will need to include:
        % - training new model
        % - score and show/plot OOB error
        
        
    end

%% NEXT BLOCK CALLBACK
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nextblock(~,~)
        % CLEAR AXES
        %         keyboard;
        try
            tmph = get(handles.f,'children');
            for uu=1:max(size(tmph))
                if strcmp(tmph(uu).Type,'axes')
                    cla(tmph(uu));
                    delete(tmph(uu));
                end
            end
            delete(handles.LFPfig);
            if handles.use_EMG
                delete(handles.EMGfig);
            end
            delete(handles.feat_axes);
            drawnow;
        catch
            keyboard;
        end
        
        % save the data
        set(handles.playbackText,'string','Saving your state data.');
        if ~isfield(handles.D,'savehere')
            if isfield(handles.D,'animal')
                basename = [handles.D.animal '_STATETIMES'];
            else
                handles.D.animal = input('Animal name ?  ','s');
                basename = [handles.D.animal '_STATETIMES'];
            end
            sDir = uigetdir(cd,'Pick saving directory.');
            matfiles_inSdir = dir([sDir filesep basename '*.mat']);
            num_to_add = numel(matfiles_inSdir) + 1;
            sFile = [basename '_' num2str(num_to_add) '.mat'];
            handles.D.savehere = [sDir filesep sFile];
            handles.D.saveDir  = sDir;
        end
        statetimes    = handles.D.AutoStates;
        save(handles.D.savehere,'statetimes');
        
        % find "real" error (how many bins were corrected)
        clear stmp stmp0 stmp_2
        stmp0 = statetimes;
        rd = handles.D.rawdat;
        stmp_2 = [zeros(size(rd,1),1) rd(:,2)];
        these_idxs = find(stmp0(:,2) >= rd(1,2) & stmp0(:,2) <= rd(end,2));
        these_idxs = unique([max(these_idxs(1)-1,1); these_idxs]);
        try
            stmp = [stmp0(these_idxs,:); [stmp0(end,1) rd(end,2)+10]];
        catch
            keyboard;
        end
        %         stmp = [stmp0(stmp0(:,2)>=rd(1,2) & stmp0(:,2)<=rd(end,2),:); [stmp0(end,1) rd(end,2)+10]];
        for aa = 1:size(rd,1)
            
            thistime = rd(aa,2);
            scount = 1;
            sdone = 0;
            while ~sdone
                try
                    if thistime >= stmp(scount,2) && thistime < stmp(scount+1,2)
                        this_state = stmp(scount,1);
                        sdone = 1;
                    else
                        scount = scount + 1;
                        sdone = 0;
                    end
                catch
                    disp('nope');
                    keyboard
                end
            end
            stmp_2(aa,1) = this_state;
        end
        
        diff_bins = rd(:,1) - stmp_2(:,1);
        perc_diff = sum(diff_bins~=0)/size(rd,1) * 100;
        
        fprintf('  *** Had to correct %.1f%% of bins. ***\n\n',perc_diff);
        handles.CHANGED_BINS(handles.D.block) = perc_diff;
        
        handles.D.upsampled_statedat = stmp_2;
        
        handles.D.no_training = 0;
        
        % Add newly scored block to training data, unless previous block
        % has just been reloaded.
        if ~handles.D.no_training && handles.master_doTrain
            if handles.D.block == 1
                try
                    handles.D.training_data = [handles.D.features{handles.D.block}, handles.D.upsampled_statedat(:,1)];
                catch
                    keyboard;
                end
            else
                handles.D.training_data = [handles.D.training_data ; ...
                    [handles.D.features{handles.D.block}, handles.D.upsampled_statedat(:,1)] ] ;
            end
            
            
            fprintf('Training data size is now %u rows by %u columns.\n',...
                size(handles.D.training_data,1),size(handles.D.training_data,2));
            
            % remove NANs from training data
            if any(any(isnan(handles.D.training_data)))
                [r_nan, ~] = find(isnan(handles.D.training_data));
                unique(r_nan);
                if numel(unique(r_nan)) == size(handles.D.training_data,1) &&...
                        handles.D.block == handles.D.X_Block
                    handles.D.training_data(isnan(handles.D.training_data)) = 0;
                else
                    handles.D.training_data(r_nan,:) = [];
                end
            end
            
            % this is the size cap for the training set
            size_cap = 10000;
            if size(handles.D.training_data,1) > size_cap
                fprintf('Too big! Clipping...\n');
                size_diff = size(handles.D.training_data,1) - size_cap;
                new_training_data = handles.D.training_data;
                new_training_data(1:size_diff,:) = [];
                handles.D.training_data = [];
                handles.D.training_data = new_training_data;
                fprintf('Training data size is now %u rows by %u columns.\n',...
                    size(handles.D.training_data,1),size(handles.D.training_data,2));
            end
            
            % Save training data (over-writing on each block)
            set(handles.playbackText,'string','Saving your feature data.');
            if ~isfield(handles.D,'savefeatshere')
                if isfield(handles.D,'animal')
                    featbasename = [handles.D.animal '_FEATURES'];
                else
                    handles.D.animal = input('Animal name ?  ','s');
                    featbasename = [handles.D.animal '_FEATURES'];
                end
                featfiles_inSdir = dir([handles.D.saveDir filesep featbasename '*.mat']);
                num_to_add = numel(featfiles_inSdir) + 1;
                sFile = [featbasename '_' num2str(num_to_add) '.mat'];
                handles.D.savefeatshere = [handles.D.saveDir filesep sFile];
            end
            features_by_block   = handles.D.features;
            training_set        = handles.D.training_data;
            changed_bins        = handles.CHANGED_BINS;
            if isfield(handles,'OOB_ERROR')
                oob_error           = handles.OOB_ERROR;
            else
                oob_error           = NaN;
            end
            
            save(handles.D.savefeatshere,'features_by_block','training_set','changed_bins','oob_error');
        else
            set(handles.playbackText,'string','Skipping feature data save.');
            fprintf(['Proceeding without adding this block''s data to training set.\n',...
                'Training data size is now %u rows by %u columns.\n'],...
                size(handles.D.training_data,1),size(handles.D.training_data,2));
        end
        
        
        try
            delete(handles.LFPfig);
            drawnow;
        end
        
        etemp = 5;
        
        while etemp == 5
            handles.D.block = handles.D.block+1;
            
            if handles.D.LFPdata(handles.D.block).duration > 300; % don't bother coding blocks of less than 5 minutes
                etemp = 55;
            end
            
            if handles.D.block == size(handles.D.LFPdata,2) && handles.D.LFPdata(handles.D.block).duration < 300
                set(handles.playbackText,'visible','on','string','last block is a miniature horse and is not worth coding. You''re done!');
                handles.littlehorse = 1;
            end
            
        end
        
        
        
        assignin('base','handles',handles);
        
        if handles.D.block<size(handles.D.LFPdata,2) && handles.littlehorse == 0
            processData;
        elseif handles.D.block>=size(handles.D.LFPdata,2) && handles.littlehorse == 0
            
            set(handles.playbackText,'visible','on','string','You made it... What does that say about you? Data are already saved.');
        end
        
    end

%% NEXT BLOCK (No training)
% this is the same as the nextblock function above, but it does not add the
% current block to the dataset. This is useful if the data in current block
% is corrupted, or if the current block includes a period where the animal
% was unplugged (e.g. for MD surgery)
    function nextblock_notrain(~,~)
        % CLEAR AXES
        %         keyboard;
        try
            tmph = get(handles.f,'children');
            for uu=1:max(size(tmph))
                if strcmp(tmph(uu).Type,'axes')
                    cla(tmph(uu));
                    delete(tmph(uu));
                end
            end
            delete(handles.LFPfig);
            if handles.use_EMG
                delete(handles.EMGfig);
            end
            delete(handles.feat_axes);
            drawnow;
        catch
            keyboard;
        end
        
        set(handles.playbackText,'string','Saving your data (without adding to training set)');
        if ~isfield(handles.D,'savehere')
            if isfield(handles.D,'animal')
                basename = [handles.D.animal '_STATETIMES'];
            else
                handles.D.animal = input('Animal name ?  ','s');
                basename = [handles.D.animal '_STATETIMES'];
            end
            sDir = uigetdir(cd,'Pick saving directory.');
            matfiles_inSdir = dir([sDir filesep basename '*.mat']);
            num_to_add = numel(matfiles_inSdir) + 1;
            sFile = [basename '_' num2str(num_to_add) '.mat'];
            handles.D.savehere = [sDir filesep sFile];
            handles.D.saveDir  = sDir;
        end
        statetimes    = handles.D.AutoStates;
        save(handles.D.savehere,'statetimes');
        
        try
            delete(handles.LFPfig);
            drawnow;
        end
        
        etemp = 5;
        
        while etemp == 5
            handles.D.block = handles.D.block+1;
            
            if handles.D.LFPdata(handles.D.block).duration > 300 % don't bother coding blocks of less than 5 minutes
                etemp = 55;
            end
            
            if handles.D.block == size(handles.D.LFPdata,2) && handles.D.LFPdata(handles.D.block).duration < 300
                set(handles.playbackText,'visible','on','string','last block is a miniature horse and is not worth coding. You''re done!');
                handles.littlehorse = 1;
            end
            
        end
        
        % if this is set to 1, the following block's data will not be added
        % to the training set. Both NEXT_BLOCK callbacks should have this
        % line re-setting no_training to 0.
        handles.D.no_training = 0;
        
        assignin('base','handles',handles);
        
        if handles.D.block<size(handles.D.LFPdata,2) && handles.littlehorse == 0
            processData;
        elseif handles.D.block>=size(handles.D.LFPdata,2) && handles.littlehorse == 0
            
            set(handles.playbackText,'visible','on','string','You made it... What does that say about you? Data are already saved.');
        end
        
    end


%% Other callback functions
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GO TO BLOCK
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function GoToBlock(~,~)
        
        e = str2num(get(handles.gotoblockBox,'string'));
        handles.D.block = e;
        
        f = [];
        g = [];
        h = [];
        
        if length(handles.D.features) < (handles.D.block) || handles.D.prior
            handles.D.no_training = 1;
        end
        
        
        assignin('base','handles',handles);
        processData;
        
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% CHANGE STATE
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeState (~,~)
        % This function allows the user to manually correct errors in the
        % automatic recognition of arousal states. The user will click on
        % the begining (x axis) and end (x axis) of the error and then
        % input the proper state code. The script below will update the
        % display and overwrite the saved data to reflect the edited
        % states.
        
        % Get user input x axis values for the onset and offset of the
        % error state code, to be replaced by input value shortly...
        
        set(handles.playbackText,'visible','on','string','Click the START of change.');
        [a1, ~] = ginput(1);
        pause (0.5);
        set(handles.playbackText,'string','Click the END of change.');
        [a2, ~] = ginput(1);
        set(handles.playbackText,'string','Enter state code. REM (1) NREM (2) Active (4) Quiet (5)');
        set(handles.changeStateinput,'enable','on');
        uicontrol(handles.changeStateinput);
        
        b = handles.D.pX - a1;
        b = abs(b);
        b = find(b == min(b));
        
        c = handles.D.pX - a2;
        c = abs(c);
        c = find(c == min(c));
        
        % Ask the user for the new state code to replace the error:
        stateoptions = [1 2 4 5];
        s = 0;
        while s == 0
            t = str2double(get(handles.changeStateinput,'string'));
            
            if any(t == stateoptions)
                s = 1;
                set([handles.changeStateinput],'Enable','off');
                set([handles.changeStateinput handles.playbackText] ,'string','');
            else
                s = 0;
                pause (0.5);
            end
        end
        
        % figure out what the next state is whose onset is determined by a2. This
        % is often a timestamp in the middle of an ongoing block, so a new state
        % will be created.
        tmp = find(handles.D.AutoStates(:,2)==a2);
        
        if isempty(tmp)
            nextstate = handles.D.AutoStates(find(handles.D.AutoStates(:,2)<a2,1,'last'),1);
        else
            nextstate = handles.D.AutoStates(tmp,1);
        end
        
        % delete the time stamps and states to be overwritten, then write
        % new entries and sort.
        handles.D.AutoStates (handles.D.AutoStates(:,2)>=a1 & handles.D.AutoStates(:,2)<=a2, : )  = [];
        
        new1 = [t a1];
        new2 = [nextstate a2];
        handles.D.AutoStates = [handles.D.AutoStates; new1; new2];
        handles.D.AutoStates = sortrows(handles.D.AutoStates,2);
        
        [~,a1_bin] = min(abs(handles.D.rawdat - a1));
        [~,a2_bin] = min(abs(handles.D.rawdat - a2));
        a2_bin = a2_bin - 1;
        
        handles.D.rawdat(a1_bin:a2_bin,1) = t;
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Plot the algorithmically coded video:
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tdat    = handles.D.AutoStates(handles.D.AutoStates(:,2)>=handles.D.pX(1) & handles.D.AutoStates(:,2)<=handles.D.pX(end),:);
        tdat = [tdat; [9 handles.D.pX(end)]];
        
        Ay = 8;
        my_ymin = -2;
        axes(handles.feat_axes);
        
        clr = [0.5 0.5 0.5; 0.2 0.2 0.3; NaN NaN NaN; 0.4 1.0 0.2; 1 1 0];
        
        for ee = 1:size(tdat,1)-1
            t0      = tdat(ee,2);
            t1      = tdat(ee+1,2);
            
            try
                yyaxis right
                rectangle( 'Position', [t0, Ay, t1 - t0 , 1],...
                    'facecolor',clr(tdat(ee,1),:),'linestyle','none');
            catch
                disp('Caught in plotting auto output');
                keyboard
            end
        end
        
        yyaxis right
        set(gca,'ylim',[my_ymin Ay+1]);
        yyaxis left
        set(gca,'ylim',[0 1]);
        set(gca,'TickLength',[0 0.025]);
        
        
        assignin('base','handles',handles);
        
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% CHECK VIDEO
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function checkvidplayback(~,~)
        % this function lets the user click on the GUI and pull up video of
        % that section of data. This is essential to differentiate between
        % active and quiet wake, and useful to correctly code
        % micro-arousals and other ambiguous sections of data
        
        file_1 = [];
        file_2 = [];
        csv_idx = [];
        csv_file = [];
        avi_file = [];
        vidread = [];
        g = [];
        h = [];
        j = [];
        k = [];
        l = [];
        m = [];
        t0 = [];
        t1 = [];
        
        % This is the new bit (1/29/2019) - atp
        % Use the 15-minute long video files instead of the one huge AVI
        % file. First find the write time of all CSV files in the relevant
        % folder. Use this to find the relevant CSV/AVI file combination
        % for video lookup.
        if ~isfield(handles.D,'CSVdir')
            csv_dir = uigetdir('C:\','Select folder containing the video and CSV files.');
            handles.D.CSVdir = csv_dir;
            handles.D.ROIcounter = 0;
            [csv_t,csv_f] = getCsvFileTimes(handles.D.CSVdir);
            handles.D.CSV_times = csv_t;
            handles.D.CSV_files = csv_f;
            assignin('base','handles',handles);
        end
        if ~isfield(handles.D,'CSV_times') || ~isfield(handles.D,'CSV_files')
            [csv_t,csv_f] = getCsvFileTimes(handles.D.CSVdir);
            handles.D.CSV_times = csv_t;
            handles.D.CSV_files = csv_f;
            assignin('base','handles',handles);
        end
        
        
        %         if ~isfield(handles.D,'AVIfile')
        %             [a,b] = uigetfile('.avi','Select the AVI file');
        %             handles.D.AVIfile = [b a];
        %             handles.D.ROIcounter = 0;
        %             handles.D.frameTimes_file = [b a(1:7) '_frameTimes.mat'];
        %             load_ft = load(handles.D.frameTimes_file);
        %             handles.D.frameTimes = load_ft.frametimes;
        %             handles.D.frameTimes = load_ft.frametimes - handles.D.hour_offset*3600;
        %             clear a b
        %             assignin('base','handles',handles);
        %         end
        
        
        %         if ~isfield(handles.D,'frameTimes')
        %             [a,b] = uigetfile('.mat','Select the frameTimes file');
        %             handles.D.frameTimes_file = [b a];
        %             load_ft = load(handles.D.frameTimes_file);
        %             handles.D.frameTimes = load_ft.frametimes - handles.D.hour_offset*3600;
        %             clear a b
        %             assignin('base','handles',handles);
        %         end
        
        
        set(handles.playbackText,'visible','on','string','Select playback START POINT');
        [t0, ~] = ginput(1);
        pause (0.5);
        set(handles.playbackText,'string','Select playback END POINT');
        [t1, ~] = ginput(1);
        
        %PRINT VERTICAL LINES ON T0 AND T1
        % dont know what happened to doing this... would be a useful
        % addition
        
        set(handles.playbackText,'string','Loading video');
        drawnow;
        
        
        % find the csv file corresponding to the chosen playback times
        file_1 = find(handles.D.CSV_times(:,1) >= t0, 1, 'first');
        file_2 = find(handles.D.CSV_times(:,1) >= t1, 1, 'first');
        
        if (file_2-file_1) > 0
            fprintf('\n\n\t*** Video spans more than one file. ***\n\n');
            %             keyboard;
            
            % VIDEO 1
            csv_idx1 = handles.D.CSV_times(file_1,2);
            csv_file1 = [handles.D.CSV_files(csv_idx1).folder filesep handles.D.CSV_files(csv_idx1).name];
            avi_file1 = [csv_file1(1:end-4) '.avi'];
            
            frametimes1 = csvread(csv_file1) - handles.D.hour_offset*3600;
            g1 = find(frametimes1 > t0,1,'first'); % first frame
            h1 = find(frametimes1 > t1,1,'first');    %last frame
            if isempty(h1)
                h1 = size(frametimes1,1);
            end
            vidread1 = VideoReader(avi_file1);
            nframes1 = round(vidread1.FrameRate * vidread1.Duration);
            if h1 > nframes1
                % this is because an error in the CSV files saving. Usually
                % what has happened is that this file has duplicate frame
                % times from the previous file/video pair. Load the
                % previous file and delete the duplicates.
%                 keyboard;
                fname_split1 = regexp(csv_file1,'_','split');
                fname_last1 = fname_split1{end};
                dot_split1 = regexp(fname_last1,'\.','split');
                file_id1 = round(str2double(dot_split1{1}))-1;
                fname_split_prev1 = fname_split1;
                fname_split_prev1{end} = [num2str(file_id1) '.csv'];
                load_prev1 = strjoin(fname_split_prev1,'_');
                if exist(load_prev1,'file')
                    frtimes_prev = csvread(load_prev1) - handles.D.hour_offset*3600;
                end                
                % removes duplicates and make a new frametimes1 array
                frtemp = frametimes1;
                frametimes1 = [];
                frametimes1 = setdiff(frtemp, frtimes_prev);
                % now recalculate the video frame boundaries
                g1 = find(frametimes1 > t0,1,'first'); % first frame
                h1 = find(frametimes1 > t1,1,'first');    %last frame
                if isempty(h1)
                    h1 = size(frametimes1,1);
                end
            end
            
            try
                vid_playback1 = read(vidread1,[g1 h1]);
            catch
                keyboard;
            end
            % VIDEO 2
            csv_idx2 = handles.D.CSV_times(file_2,2);
            csv_file2 = [handles.D.CSV_files(csv_idx2).folder filesep handles.D.CSV_files(csv_idx2).name];
            avi_file2 = [csv_file2(1:end-4) '.avi'];
            
            frametimes2 = csvread(csv_file2) - handles.D.hour_offset*3600;
            g2 = find(frametimes2 > t0,1,'first'); % first frame
            h2 = find(frametimes2 > t1,1,'first');    %last frame
            vidread2 = VideoReader(avi_file2);
            
            nframes2 = round(vidread2.FrameRate * vidread2.Duration);
            if h2 > nframes2
                % this is because an error in the CSV files saving. Usually
                % what has happened is that this file has duplicate frame
                % times from the previous file/video pair. Load the
                % previous file and delete the duplicates.
%                 keyboard;
                fname_split2 = regexp(csv_file2,'_','split');
                fname_last2 = fname_split2{end};
                dot_split2 = regexp(fname_last2,'\.','split');
                file_id = round(str2double(dot_split2{1}))-1;
                fname_split_prev2 = fname_split2;
                fname_split_prev2{end} = [num2str(file_id) '.csv'];
                load_prev2 = strjoin(fname_split_prev2,'_');
                if exist(load_prev2,'file')
                    frtimes_prev2 = csvread(load_prev2) - handles.D.hour_offset*3600;
                end
                % removes duplicates and make a new frametimes1 array
                frtemp2 = frametimes2;
                frametimes2 = [];
                frametimes2 = setdiff(frtemp2, frtimes_prev2);
                % now recalculate the video frame boundaries
                g2 = find(frametimes2 > t0,1,'first'); % first frame
                h2 = find(frametimes2 > t1,1,'first');    %last frame
            end
            
            try
                vid_playback2 = read(vidread2,[g2 h2]);
            catch
                keyboard;
            end
            % play them back
            close(figure(100));
            close(figure(200));
            
            if isempty(handles.D.ROI)
                %             first_frame = read(vidread, g);
                %             figure(200); imshow(g);
                %             tmp_rect = imrect();
                %             handles.D.ROI = tmp_rect;
                handles.D.ROIcounter = 1000;
            end
            
            if handles.D.ROIcounter < 10
                h = figure(100);
                inew = vid_playback(:,:,:,1).*uint8(repmat(handles.D.ROI,[1,1,3]));
                imshow(inew)
                handles.D.ROIcounter = handles.D.ROIcounter + 1;
            end
            %         imshow(f(:,:,:,1));
            
            k = figure('NumberTitle','off','Name','CLOSE THIS WHEN DONE WITH VIDEO 1');
            implay(vid_playback1,30);
            set(findall(0,'tag','spcui_scope_framework'),'units','normalized','position',[0.3 0.3 0.4 0.4]);
            uiwait(k);
            
            implay(vid_playback2,30);
            set(findall(0,'tag','spcui_scope_framework'),'units','normalized','position',[0.3 0.3 0.4 0.4]);
            
            
        else
            fprintf('Found video file.\n');
            
            
            csv_idx = handles.D.CSV_times(file_1,2);
            csv_file = [handles.D.CSV_files(csv_idx).folder filesep handles.D.CSV_files(csv_idx).name];
            avi_file = [csv_file(1:end-4) '.avi'];
            fprintf('Opening video file: %s\n',avi_file);
            
            frametimes = csvread(csv_file) - handles.D.hour_offset*3600;
            g = find(frametimes > t0,1,'first'); % first frame
            h = find(frametimes > t1,1,'first');    %last frame
            vidread = VideoReader(avi_file);
            nframes = vidread.Duration * vidread.FrameRate;
            
            % this can happen if there is an error in the video acquisition
            % resulting in overlaps in frametimes in the csv files
            if h > nframes
                % load previous csv files to check for overlap
                % this is because an error in the CSV files saving. Usually
                % what has happened is that this file has duplicate frame
                % times from the previous file/video pair. Load the
                % previous file and delete the duplicates.
%                 keyboard;
                fname_split = regexp(csv_file,'_','split');
                fname_last = fname_split{end};
                dot_split = regexp(fname_last,'\.','split');
                file_id = round(str2double(dot_split{1}))-1;
                fname_split_prev = fname_split;
                fname_split_prev{end} = [num2str(file_id) '.csv'];
                prev_csv = strjoin(fname_split_prev,'_');
                
%                 prev_csv = [handles.D.CSV_files(csv_idx-1).folder filesep handles.D.CSV_files(csv_idx-1).name];
                prev_frametimes = csvread(prev_csv) - handles.D.hour_offset*3600;
                [~,bad_idx,~] = intersect(frametimes,prev_frametimes);
                frametimes(bad_idx) = [];
                g = find(frametimes > t0,1,'first'); % first frame
                h = find(frametimes > t1,1,'first');    %last frame
            end
            
            
            try
                vid_playback = read(vidread,[g h]);
            catch
                fprintf("***Could not open video***\n");
                keyboard;
            end
            
            
            %         addpath(b);
            %         c = VideoReader(handles.D.AVIfile);
            %         d = find(handles.D.frameTimes>t0,1,'first'); % first frame
            %         e = find(handles.D.frameTimes>t1,1,'first');    %last frame
            %         f = read(c,[d e]);
            
            close(figure(100));
            close(figure(200));
            
            if isempty(handles.D.ROI)
                %             first_frame = read(vidread, g);
                %             figure(200); imshow(g);
                %             tmp_rect = imrect();
                %             handles.D.ROI = tmp_rect;
                handles.D.ROIcounter = 1000;
            end
            
            if handles.D.ROIcounter < 10
                h = figure(100);
                inew = vid_playback(:,:,:,1).*uint8(repmat(handles.D.ROI,[1,1,3]));
                imshow(inew)
                handles.D.ROIcounter = handles.D.ROIcounter + 1;
            end
            %         imshow(f(:,:,:,1));
            
            implay(vid_playback,30);
            set(findall(0,'tag','spcui_scope_framework'),'units','normalized','position',[0.3 0.3 0.4 0.4]);
        end
        set(handles.playbackText,'string',[]);
    end

% the following functions are useful if you want to stop training for a few
% blocks (e.g. if the model is coding very well you might want to keep it
% as is), resume the training, or change the ML algorithm used for
% classification (NOTE: only RF is implemented as of 6/26/2020).
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% STOP TRAINING
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function stopTraining(~,~)
        handles.master_doTrain = 0;
        set(handles.trainBox,'string','No');
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% RESUME TRAINING
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function resumeTraining(~,~)
        handles.master_doTrain = 1;
        set(handles.trainBox,'string','Yes');
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% CHANGE ALGORITHM
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function changeAlgorithm(~,~)
        
        handles.D.ml_algorithm = [];
        new_algo = handles.mlAlgoBox.String{handles.mlAlgoBox.Value};
        handles.D.ml_algorithm = new_algo;
        switch new_algo
            case 'RF'
                algo_text = 'Random Forest';
            case 'ECOC'
                algo_text = 'Error-correcting output code';
        end
        set(handles.playbackText,'string',sprintf('Using %s algorithm.',algo_text));
        
    end

end