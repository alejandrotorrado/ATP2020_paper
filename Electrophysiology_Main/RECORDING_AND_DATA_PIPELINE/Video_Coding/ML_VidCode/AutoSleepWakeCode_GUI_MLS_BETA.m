%% BIG NOTE
% RF may not be the best choice for this problem.
% Maybe look into Error Correcting Output Codes (ecoc)

function AutoSleepWakeCode_GUI_MLS_BETA
warning('off');

disp('Initializing GUI');

set(0,'Units','normalized')

%% GUI elements
% muscle beach logo
% mb_logo = imread('/Volumes/turrigiano-lab/ATP_MAIN/misc_stuff/musclebeach_logo.tiff');

% figure color
f_color = [189, 243, 250]./255;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Main figure and loading buttons

handles.f = figure('WindowStyle','normal','HitTest','on','Visible','on','NumberTitle',...
    'off','Color',f_color,'Units','normalized','Position',...
    (get(0,'ScreenSize')-[-0.07 -0.07 0.1 0.1]));

handles.priortext = uicontrol('Style','text','visible','on','Backgroundcolor',f_color,...
    'String','Load prior data?','units','normalized','Position',[0.42 0.83 0.1 0.04],...
    'Horizontalalignment','center','FontSize',16,'fontname','Times');

handles.loadprioryes = uicontrol('Style','pushbutton','String','Yes',...
    'units','normalized','Position', [0.41 0.78 0.05 0.04],'Callback',@LoadPrior,...
    'fontname','Times','FontSize',15,'Visible','on');

handles.loadpriorno = uicontrol('Style','pushbutton','String','No',...
    'units','normalized','Position', [0.49 0.78 0.05 0.04],'Callback',@LoadPrior,...
    'fontname','Times','FontSize',15,'Visible','on');

handles.loadLFP_button = uicontrol('Style','pushbutton','String','Load LFP data',...
    'units','normalized','Position', [0.42 0.58 0.2 0.1],'Callback',@LoadLFPdata,...
    'fontname','Times','FontSize',16,'Visible','off');

handles.EMGdir_button = uicontrol('Style','pushbutton','String','Select EMG directory',...
    'units','normalized','Position', [0.42 0.58 0.2 0.1],'Callback',@EMGdir,...
    'fontname','Times','FontSize',16,'Visible','off');

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

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% playback text

% logo_axes = axes('position',[.9 .01 .085 .085]);
% logo_image = image(mb_logo(:,:,1:3));
% logo_image.AlphaData = 0.4;
% set(logo_axes,'box','off','xtick',[],'xticklabel',[],...
%     'ytick',[],'yticklabel',[],'Ycolor',f_color,'Xcolor',f_color,...
%     'color',f_color);

%% CALLBACK FUNCTION DEFINITIONS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function LoadPrior(~,eventdata)
     % If the user has already done some coding, load the previous work and
     % make a note that this exists so as to (later) start at the correct
     % point in time.
       
       if strcmp(eventdata.Source.String,'Yes')
           
           [sFile, sDir] = uigetfile('.mat','Pick your previous statetimes.');
           
           % find number of file
           underscore_split = regexp(sFile,'_','split');
           handles.D.animal = underscore_split{1};
           matsplit = regexp(underscore_split{end},'\.','split');
           file_ID = matsplit{1};
           
           fFile = [handles.D.animal '_FEATURES_' file_ID '.mat'];
           
           statedat = load([sDir sFile]);
           featdat = load([sDir fFile]);
           
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
       
       set(handles.loadLFP_button,'visible','on');
       set([handles.loadprioryes handles.loadpriorno handles.priortext],'visible','off');
       assignin('base','handles',handles)
        
    end


    function LoadLFPdata(~,~)
        
        set(handles.loadLFP_button,'visible','off');
        [eFile, eDir] = uigetfile('','Select the LFP data file');
        set(handles.playbackText,'String','Loading LFP data. This will take a few minutes.','visible','on');
        drawnow
        temp = load([eDir eFile]);
        set(handles.playbackText,'String','','visible','off');
        set(handles.EMGdir_button,'visible','on');
        
        handles.D.LFPdata = temp.LFPinfo;
        LFPstarts = [handles.D.LFPdata.startTime];
        
        
        % if LFP start time is NOT GMT corrected - do this.
        %         if ~isfield(handles.D.LFPdata,'GMTcor');
        %             for rr = 1:size(handles.D.LFPdata,2);
        %                 handles.D.LFPdata(rr).startTime    =  handles.D.LFPdata(rr).startTime - 3600*4;
        %                 handles.D.LFPdata(rr).GMTcor       = 1;
        %             end
        %         end
        
        if handles.D.prior == 0
            handles.D.block = 1; % CHANGE BACK TO 1
        else

            a = [handles.D.LFPdata.startTime]';
            b = find(a>handles.D.statetimes(end,2),1);
            handles.D.block = b;

        end
        handles.D.AutoStates = [];
        handles.littlehorse = 0;
        handles.master_doTrain = 1;
        
        assignin('base','handles',handles);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
    end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function EMGdir(~,~)
        
        set(handles.EMGdir_button,'visible','off');
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
        
        
        EMGindex = zeros(size(handles.D.enames,2),2);
        for ee = 1:size(handles.D.enames,2);
            
            t1 = []; t2 = [];
            
            EMGindex(ee,1) = ee;
            t1 = strfind(handles.D.enames{ee},'_EMG');
            t2 = strfind(handles.D.enames{ee},'.mat');
            EMGindex(ee,2) = str2double(handles.D.enames{ee}(t1+4:t2-1));
            
        end
        
        EMGindex = sortrows(EMGindex,2);
        
        handles.D.EMGidx = EMGindex;
        
        set(handles.loadMvmt_button,'visible','on');
        assignin('base','handles',handles);
    end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function LoadMvmt(~,~)
        
        %  Load automated movement analysis - - - - - - - - - - - - - - - -
        
        [eFile, eDir] = uigetfile('','Select movement data');
        
        temp = load([eDir eFile]);
        if isfield(temp,'DATA')
            
            pmove = temp.DATA.smooth_movement';
            %             pframes = temp.DATA.frame_times;
            %             if size(pmove,1) ~= size(pframes,1)
            %                 szdiff = abs(size(pmove,1) - size(pframes,1));
            %                 if szdiff == 1
            %                     if size(pframes,1) > size(pmove,1)
            %                         pframes(1) = [];
            %                     else
            %                         pmove(1) = [];
            %                     end
            %                 else
            %                     error('Number of frames and number of frame timestamps are not matching!');
            %                 end
            %             end
            try
                handles.D.movement  = [pmove temp.DATA.frame_times(1:end)];
            catch
                keyboard
            end
            handles.D.ROI       = temp.DATA.mask;
            
            % check for GMT time.
            lfp_start = handles.D.LFPdata(1).startTime;
            lfp_end = handles.D.LFPdata(end).startTime;
            mvt_start = handles.D.movement(1,2);
            mvt_end = handles.D.movement(end,2);
            
            
            
            if abs(lfp_start-mvt_start) >= (4*3600 - 10) && ...
                    abs(lfp_end-mvt_end) >= (4*3600 - 10)
                
                warn_Str = sprintf(['WARNING! It looks like the LFP times and movement data',...
                    ' times are misaligned. Possibly due to missing GMT time correction.\n',...
                    'The start times are:\n%s for LFP\n%s for movement.\n\n',...
                    'The end times are:\n%s for LFP\n%s for movement.\n\n',...
                    '\n_____     *** Fixing GMT offset ***     _____\n\n\n'],...
                    datestr(unixtime(lfp_start)),datestr(unixtime(mvt_start)),...
                    datestr(unixtime(lfp_end)),datestr(unixtime(mvt_end)));
                set(handles.GMTplaybackText,'visible','on','string',warn_Str);
                pause(2.5);
                
                % Movement timestamps are in GMT time. Adjust this to match
                % LFP timestamps.
                first_timestamp = unixtime(handles.D.movement(1,2));
                in_local_time = TimezoneConvert_ATP(first_timestamp,'UTC','America/New_York');
                in_local_unixtime = unixtime(datevec(in_local_time)); 
                hour_offset = round(abs(handles.D.movement(1,2) - in_local_unixtime))/3600;
                handles.D.movement(:,2) = handles.D.movement(:,2) - hour_offset*3600;
                
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
            try
                handles.D.ROI = temp.outdata.ROI;
            catch
                handles.D.ROI = temp.outdata.DATA.mask;
            end
        end
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        set(handles.loadMvmt_button,'visible','off');
        assignin('base','handles',handles);
        processData;
    end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%                         Finished UI Load Data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function processData(~,~)
        
        set(handles.playbackText,'string','Processing data.');
        drawnow;
        try
            cla(handles.LFPfig);
        end
        
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
        
        % New method:
        % Get absolute power from Power Spectral Density, for delta,
        % theta, and whole range.
        delta_absolute_power = sum(handles.D.P(handles.D.F >= delta_freq(1) & handles.D.F <= delta_freq(2),:));
        theta_absolute_power = sum(handles.D.P(handles.D.F >= theta_freq(1) & handles.D.F <= theta_freq(2),:));
        total_absolute_power = sum(handles.D.P(:,:));
        
        % Get relative power fraction for each frequency band and smooth
        delta_relative_power = delta_absolute_power ./ total_absolute_power;
        deltapower = smooth(delta_relative_power,200);
        delta_mean = nanmean(deltapower);
        delta_std  = nanstd(deltapower);
        delta_zscore = (deltapower - delta_mean) ./ delta_std;
        
        theta_relative_power = theta_absolute_power ./ total_absolute_power;
        thetapower = smooth(theta_relative_power,200);
        theta_mean = nanmean(thetapower);
        theta_std  = nanstd(thetapower);
        theta_zscore = (thetapower - theta_mean) ./ theta_std;
        
        handles.D.deltapower = deltapower;
        handles.D.thetapower = thetapower;
        handles.D.delta_z    = delta_zscore;
        handles.D.theta_z    = theta_zscore;
        
        
        %% ATP NOTE 01
        % Load and process EMG signal.
        
        % Load EMG data
        etemp = load([handles.D.eDir  filesep handles.D.enames{ handles.D.EMGidx(handles.D.block,1) } ]);
        emgLtemp = etemp.EMGdata.f_EMG_L;
        emgRtemp = etemp.EMGdata.f_EMG_R;
        
        % old method:
        %{
        if var(emgLtemp) > var(emgRtemp);
            handles.D.emgData = emgLtemp;
        else
            handles.D.emgData = emgRtemp;
        end
        %}
        
        % Average together the signal from both hemispheres.
        emgavg = nanmean([emgLtemp; emgRtemp]);
        handles.D.emgData = emgavg;
        
        % Smooth and resample EMG signal
        emgX = linspace(handles.D.LFPdata(handles.D.block).startTime,  handles.D.LFPdata(handles.D.block).startTime...
            + handles.D.LFPdata(handles.D.block).duration, size(handles.D.emgData,2) );
        emgY = handles.D.emgData;
        
        smooth_span = 2000;
        handles.D.semgY = smooth(emgY,smooth_span);
        
        % take a moving average to resample the EMG data (to match LFP length)
        bb = round(linspace(1, size(handles.D.semgY,1), size(handles.D.pX,2)+1  ));
        
        for ee = 1:size(bb,2)-1
            handles.D.rs_emg(ee) = nanmean( handles.D.semgY( bb(ee):bb(ee+1) , 1 ));
        end
        rs_emg = handles.D.rs_emg;
        
        % Transform EMG to z-score for generalizability. NOTE: this is
        % z-score for this block. May need to extend to whole experiment to
        % make truly generalizable.
        emg_mean = nanmean(rs_emg);
        emg_std  = nanstd(rs_emg);
        emg_zscore = (rs_emg - emg_mean) ./ emg_std;
        handles.D.emg_z = emg_zscore;
        
        
        
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

        % remove unwanted times
        prev_times = find(mvmt(:,2) < handles.D.pX(1));
        post_times = find(mvmt(:,2) > handles.D.pX(end));
        mvt_zscore([prev_times; post_times],:)   = [];
%         mvt_zscore(mvmt(:,2)>handles.D.pX(end),:) = [];
        
        bb = round(linspace(1, size(mvt_zscore,1), size(handles.D.pX,2)+1  ));
        
        for ee = 1:size(bb,2)-1
            rs_mvt_z(ee) = nanmean( mvt_zscore( bb(ee):bb(ee+1) , 1));
        end
        
        % adjust so that no movement = 0
        rs_mvt_z = rs_mvt_z + abs(min(rs_mvt_z));
        handles.D.rs_mvt_z = rs_mvt_z;
        
        %% ATP NOTE 06
        % FEATURES TO BE USED
        % _ EMG
        % _ Delta
        % _ Theta
        % _ Mvt
        % _ (Delta - Theta)
        % _ var(EMG)
        % _ var(Mvt)
        
         % Theta-Delta difference
        dt_diff = deltapower - thetapower;
        handles.D.dt_diff = dt_diff;
        
        if size(rs_mvt_z,1) ~= size(deltapower,1) && size(rs_mvt_z,1) == size(deltapower,2)
            rs_mvt_z = rs_mvt_z';
        end
        if size(thetapower,1) ~= size(deltapower,1) && size(thetapower,1) == size(deltapower,2)
            thetapower = thetapower';
        end
        if size(emg_zscore,1) ~= size(deltapower,1) && size(emg_zscore,1) == size(deltapower,2)
            emg_zscore = emg_zscore';
        end
        if size(dt_diff,1) ~= size(deltapower,1) && size(dt_diff,1) == size(deltapower,2)
            dt_diff = dt_diff';
        end
        
        % COMPILE DATA REQUIRED FOR FEATURE CALCULATION IN A MATRIX
        handles.D.feature_data = [deltapower, thetapower, dt_diff, rs_mvt_z, emg_zscore];
             
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
        
        testdat = [];
        if handles.D.block <= handles.D.X_Block
            sampling_rate = 2; % Hz
            set(handles.playbackText,'visible','on','string','Scoring data');
            [testdat, rawdat, block_features] = tensectest_new(handles.D.feature_data, handles.D.pX, sampling_rate);
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
            %  - 'ECOC' : Error-correcting output code
            
            if handles.master_doTrain
                handles.D.trained_Mdl = [];
                
                % train model based on previous blocks
                
                set(handles.playbackText,'string',sprintf('Training and saving %s model...',handles.D.ml_algorithm));
                drawnow;
                talg0 = tic;
                [trained_mdl,oob_err,mdl_path] = AutoVidCode_trainModel(handles.D.training_data,...
                    handles.D.animal,handles.D.ml_algorithm,[],0);
                set(handles.oobBox,'String',sprintf('%.2f%%',100*oob_err));
                
                %--------------------------------------------------------------
                % ******________________ IMPORTANT NOTE: ________________******
                % By default the above function uses n=200 trees and does not
                % display anything. To change the number of trees, pass it as
                % the 4th argument to the function. To display the OOB error as
                % a function of n_Trees, pass display_flag = 1 as the 5th
                % argument.
                %--------------------------------------------------------------
                %--------------------------------------------------------------
                
                % THIS SHOULD ALWAYS BE 1!!! 0 ONLY FOR DEBUGGING
                do_save_mdl = 1;
                if do_save_mdl
                    Mdl_struct.Mdl = trained_mdl;
                    Mdl_struct.ML_algo = handles.D.ml_algorithm;
                    Mdl_struct.block = handles.D.block;
                    
                    % Saving model
                    fprintf('\n  Saving your ML model...\n');
                    ts0 = tic;
                    save(mdl_path,'Mdl_struct');
                    ts1 = toc(ts0);
                    fprintf('Done! Time elapsed: %.2f seconds.\n\n');
                end
                
                talg1 = toc(talg0);
                set(handles.playbackText,'string',sprintf('Training and saving %s model... That took %.2f seconds.',...
                    handles.D.ml_algorithm,talg1));
                handles.D.trained_Mdl = trained_mdl;
            else
                set(handles.playbackText,'string',sprintf('Training is OFF! Using last trained model.'));
            end
            pause(.1);
            
%             loaded_Mdl = AutoVidCode_loadModel(handles.D.animal,ml_algorithm);
            
            % parameters for video scoring. Bin size for scoring and
            % sampling rate of input data
            bin_sz = 10; % seconds
            sampling_rate = 2; % Hz
            set(handles.playbackText,'visible','on','string',...
                sprintf('Auto-scoring using %s model. Please wait...',handles.D.ml_algorithm));
            % Call to autoscoring function
            [testdat, rawdat, block_features] = auto_video_score(handles.D.feature_data, handles.D.trained_Mdl, ...
                handles.D.pX, bin_sz, sampling_rate);
            handles.D.features{handles.D.block} = block_features;
        end
        
        handles.D.rawdat = rawdat;
           
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
        
        %% ATP NOTE 09
        % Review the post-processing below.
        
        % find instances of REM preceeded by anything other than NREM and
        % fix if appropriate:
        a = testdat(3:end-2,1); b = testdat(2:end-3,1); c = testdat(1:end-4,1);
        tmp = find(a == 1 & b > 3 & c > 3); % this will return indices of rem preceeded by 20sec of quiet waking
        
        
        try
        for ee = 1:size(tmp,1)
            
            if tmp(ee) == 1 && testdat(1) > 3
                
                testdat(1:3) = 1;
                
            elseif size(a,1)>=tmp(ee)+4 && c(tmp(ee)-1) == 2 || ...
                    size(a,1)>=tmp(ee)+4 && c(max(tmp(ee)-2,1)) == 2 % if bin preceeding the 20sec of Q is NREM, convert the 20s to REM
                
                testdat(tmp(ee)-1:tmp(ee)+1,1) = 2; % remember that 'a' is offset by 2, and you're correcting the two prior bins.
                
            elseif size(a,1)>=tmp(ee)+4 && any(c(tmp(ee) - 1) == [5 4]) && any(a(tmp(ee) + 1) == [5 4]) || size(a,1)>=tmp(ee)+4 && ...
                    any(c(tmp(ee) - 1) == [5 4]) && any(a(tmp(ee) + 4) == [5 4])
                
                
                flanks = a(tmp(ee)-1 : tmp(ee) + 5);
                flanks = mode(flanks(flanks~= 1));
                
                testdat(tmp(ee)+1 : tmp(ee) + 7,1) = flanks; % correct for random 10-20sec REM bin in the middle of waking.

                 
            end
            
        end
        catch
            keyboard
        end
        
        
        %% 
        % - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - -
        %  - - - - - - - - - - - DO THE PLOTTING - - - - - - - - - - - - -
        % - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - -
        
        % Plot the spectrogram:
        try
            delete(handles.LFPfig);
            delete(handles.feat_axes);
        end
        handles.LFPfig = axes('Position',[.05 .13 .9 .75],'visible','off','Parent',handles.f);
        colormap parula
        
        axes(handles.LFPfig);
        psd_axes_pos = handles.LFPfig.Position;
        imagesc(handles.D.pX,( handles.D.F ),10*log10(abs(handles.D.P)));
        set(handles.LFPfig,'ydir','normal','visible','off')
        %xtemp = (handles.D.pX - (handles.D.pX(1)))/60;
        
        handles.feat_axes = axes('Position',psd_axes_pos,'Color','none');
        
        % old:
        %{
        plot(handles.D.pX, thetapower,'k','linewidth',2);
        plot(handles.D.pX, sm_delta,'r','linewidth',2);
        plot(mvmt(:,2),mvmt(:,1),'m--','linewidth',2);
        plot(handles.D.pX, handles.D.rs_emg,'.','linewidth',0.5);
        %}
        
        % plot the z-score of movement and EMG on right y axis
        yyaxis right; hold on;
        plot_mvt = plot(handles.D.pX, handles.D.rs_mvt_z,'m--','linewidth',2);
        plot_emg = plot(handles.D.pX, handles.D.emg_z,'.','linewidth',0.5,...
            'color',[0 .45 .74]);
        hold off;
        % plot fraction of delta and theta power on left y axis (and
        % delta-theta)
        yyaxis left; hold on;
        set(gca,'visible','on','fontsize',14);
        plot_theta = plot(handles.D.pX, handles.D.thetapower,'k','linewidth',2);
        plot_delta = plot(handles.D.pX, handles.D.deltapower,'r','linewidth',2);
        plot_dt    = plot(handles.D.pX, handles.D.dt_diff,'--c','linewidth',1.5);
        hold off;
        
        % legend
        plot_legend = legend([plot_theta,plot_delta,plot_mvt,plot_emg],...
            {'Theta','Delta','Movement','EMG'},'fontsize',14);
        set(plot_legend,'position',[0.85 0.66 0.08 0.07],'color',[1 1 1]);
        
        % set the x axis labels to minutes elapsed in the plotted dataframe
        xlims        = get(handles.feat_axes,'xlim');
        xlims2 = get(handles.LFPfig,'xlim');
        onesec_inpts = size(handles.D.pX,2)/3600;
        set(handles.feat_axes,'xlim',xlims2);
        set(handles.LFPfig,'xlim',xlims2);
        set(handles.feat_axes,'XTick',[xlims2(1):600:xlims2(2)]);
        currlabel   = get(handles.feat_axes,'XTick');
        newlabel    = currlabel - xlims2(1);
        set(gca,'XTickLabel',newlabel/60);
        
        xl          = xlabel('Time (minutes)');
        yyaxis left
        yl_1        = ylabel('Fraction of PSD (Delta and Theta features)');
        yyaxis right
        yl_2        = ylabel('Z-score (EMG and Mvt features');
        set([xl yl_1 yl_2],'Fontname','Myriad','Fontsize',16)
        
        
        
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Plot the algorithmically coded video:
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Ay = 8;
        my_ymin = -2;
        axes(handles.feat_axes);
        
        clr = [0.5 0.5 0.5; 0.2 0.2 0.3; NaN NaN NaN; 0.4 1.0 0.2; 1 1 0];
        
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
        clear stmp stmp_2
        stmp0 = statetimes;
        rd = handles.D.rawdat;
        stmp_2 = [zeros(size(rd,1),1) rd(:,2)];
        stmp = [stmp0(stmp0(:,2)>=rd(1,2),:); [stmp0(end,1) rd(end,2)+1]];
        for aa = 1:size(rd,1)
            
            thistime = rd(aa,2);
            scount = 1;
            sdone = 0;
            
            try
            while ~sdone
                if thistime >= stmp(scount,2) && thistime < stmp(scount+1,2)
                    
                    this_state = stmp(scount,1);
                    sdone = 1;
                else
                    scount = scount + 1;
                    sdone = 0;
                end
            end
            catch
                keyboard
            end
            stmp_2(aa,1) = this_state;
        end

        diff_bins = rd(:,1) - stmp_2(:,1);
        perc_diff = sum(diff_bins~=0)/size(rd,1) * 100;
        fprintf('  *** Had to correct %.1f%% of bins. ***\n\n',perc_diff);
        
        handles.D.upsampled_statedat = stmp_2;
        
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
            
            save(handles.D.savefeatshere,'features_by_block','training_set');
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
        
        handles.D.no_training = 0;
        
        assignin('base','handles',handles);
        
        if handles.D.block<size(handles.D.LFPdata,2) && handles.littlehorse == 0
            processData;
        elseif handles.D.block>=size(handles.D.LFPdata,2) && handles.littlehorse == 0
            
            set(handles.playbackText,'visible','on','string','You made it... What does that say about you? Data are already saved.');
        end
        
    end

%% NEXT BLOCK (No training)
function nextblock_notrain(~,~)
       
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
        
        a = [];
        b = [];
        c = [];
        d = [];
        e = [];
        f = [];
        g = [];
        h = [];
        t0 = [];
        t1 = [];
        
        if ~isfield(handles.D,'AVIfile');
            [a,b] = uigetfile('.avi','Select the AVI file');
            handles.D.AVIfile = [b a];
            handles.D.ROIcounter = 0;
            clear a b
            assignin('base','handles',handles);
        end
        
        
        set(handles.playbackText,'visible','on','string','Select playback START POINT');
        [t0, ~] = ginput(1);
        pause (0.5);
        set(handles.playbackText,'string','Select playback END POINT');
        [t1, ~] = ginput(1);
        
        %PRINT VERTICAL LINES ON T0 AND T1
        
        set(handles.playbackText,'string','Loading video');
        
        c = VideoReader(handles.D.AVIfile);
        d = find(handles.D.movement(:,2)>t0,1,'first'); % first frame
        e = find(handles.D.movement(:,2)>t1,1,'first');    %last frame
        f = read(c,[d e]);   
        
        close(figure(100));
        
        if handles.D.ROIcounter < 10
            h = figure(100);
            inew = f(:,:,:,1).*uint8(repmat(handles.D.ROI,[1,1,3]));
            imshow(inew)
            handles.D.ROIcounter = handles.D.ROIcounter + 1;
        end
        %         imshow(f(:,:,:,1));
        
        implay(f,30);
        set(findall(0,'tag','spcui_scope_framework'),'units','normalized','position',[0.3 0.3 0.4 0.4]);
        
        set(handles.playbackText,'string',[]);
        
    end

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