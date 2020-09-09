function timeVid_v4_2animals(obj, event, string_arg,save_int) % obj, event

%%%% ATP 11/23/2015
% Change this code:
%  1) Remove the LED tracking part
%  2) Increase the frame rate (how often a frame is acquired/saved)
%      To do this, I will remove count_master and just let the function
%      save a frame everytime it is called. This should make the frame rate
%      a lot faster.
%  3) Keep the timing pulses
%
% -------------------------------------------------------------------------
% First pass:
% Comment out all the unnecessary code - then test it.
% If it works, can delete commented out parts.
% -------------------------------------------------------------------------

global handles;
verbose = 0;

if verbose ==1, fprintf('Beginning timer function.\n'), end

handles.count_master    = handles.count_master+1;
handles.count_A         = handles.count_A+1;
% all these counters are now useless
%{
handles.smallcount_A    = handles.smallcount_A+1;
handles.count_B         = handles.count_B+1;
handles.smallcount_B    = handles.smallcount_B+1;
%}

if handles.count_master==1; %%%%% WHAT DO I DO WITH THIS?
    tic;
end
% DO NOT UNDERSTAND THIS - if (in previous version) count was more than 1,
% when would the tic begin?
% ANS: toc will just mark the elapsed time since the last tic. Which
% presumably would be a previous run of the script? Dont know if that
% works. Need to check this.
% *** For now, can set as count_master ***

% - - - - - - - - -
% this is code to produce the shared "check file" on the Rpi. currently,
% raspberry pi will check  the time stamps on this file and will get pissed
% off if they stop increasing - KBH 1/11/16
if unixtime(clock) - handles.timezero>50; % only do this after ~50s of data
    uploadchecker_pi(unixtime(clock),handles.checkerDir,0);
    handles.timezero = unixtime(clock);
end



% - - - - - - - - - - - -




b=toc;
num=obj.TasksExecuted;
interval = save_int - 0.2; % interval between save events (shifted by 200 msec)

% save data and send pulse every 1800 sec approx (1/2 hour) - make 1799.8
if b >= interval %1799.8;
    tic % this is the tic that determines the 1/2 hour save period after count_master > 1
    
    %     t0 = tic; % this is my testing tic. ATP 11/23/15
    if verbose == 1, fprintf('\nBeginning save routine.\n'), end
    
    %     all_yStart = min(handles.yStart_A,handles.yStart_B);
    %     all_yEnd = max(handles.yEnd_A,handles.yEnd_B);
    %     all_xStart = min(handles.xStart_A,handles.xStart_B);
    %     all_xEnd = max(handles.xEnd_A,handles.xEnd_B);
    
    trigger(handles.vid)
    
    datT            = clock;
    svdata           = getdata(handles.vid,1);
    
    % shouldn't need all this:
    %{
    data            = (svdata(:,:,3));
    %data            = rgb2gray(data);
    data_all        = data(all_yStart:all_yEnd,all_xStart:all_xEnd);
    data_A          = data(handles.yStart_A:handles.yEnd_A,handles.xStart_A:handles.xEnd_A);
    data_B          = data(handles.yStart_B:handles.yEnd_B,handles.xStart_B:handles.xEnd_B);
    data_all        = imresize(data_all,0.2,'bilinear');
    data_A          = imresize(data_A,0.2,'bilinear');
    data_B          = imresize(data_B,0.2,'bilinear');
    %}
    
    % don't need these filters...
    %{
    if ~isfield(handles,'filters_A') || ~isfield(handles,'filters_B');
        warning(['timeVid_v3_2animals.m :: Trying to save data but IR '...
            'filters are not defined. Calling keyboard stop (line 64).']);
        keyboard;
    end
    %}
    
    % should be able to remove alpha_v3 (and all the code that follows)
    % completely
    %    [x_A,y_A,timestamp_A] = alpha_v3 (data_A,handles.filters_A);
    %    [x_B,y_B,timestamp_B] = alpha_v3 (data_B,handles.filters_B);
    %{
    % ANIMAL A FIX 1
    if length(x_A)~=1 || length (y_A)~=1;
        disp('XY problem detected and trying to fix');
        
        if length(x_A)>1;
           x_A = x_A(1);
        end
        if length(y_A)>1;
                   y_A = y_A(1);
        end
        if isempty(x_A);
            x_A = 1;
        end
        if isempty(y_A);
            y_A = 1;
        end
    end
    
    % ANIMAL B FIX 1
    if length(x_B)~=1 || length (y_B)~=1;
        disp('XY problem detected and trying to fix');
        
        if length(x_B)>1;
           x_B = x_B(1);
        end
        if length(y_B)>1;
                   y_B = y_B(1);
        end
        if isempty(x_B);
            x_B = 1;
        end
        if isempty(y_B);
            y_B = 1;
        end
    end
    
    % ANIMAL A FIX 2
    %x = (handles.xEnd - handles.xStart) - x;
    if x_A==1 && y_A==1;
        try
            handles.XY_list_A(handles.smallcount_A,1)= handles.XY_list_A(handles.smallcount_A-1,1);
            handles.XY_list_A(handles.smallcount_A,2)= handles.XY_list_A(handles.smallcount_A-1,2);
        catch
            handles.XY_list_A(handles.smallcount_A,1)= 1;
            handles.XY_list_A(handles.smallcount_A,2)= 1;
        end
        handles.XY_list_A(handles.smallcount_A,3)= timestamp_A;
    else
        handles.XY_list_A(handles.smallcount_A,1)=[x_A];
        handles.XY_list_A(handles.smallcount_A,2)=[handles.yEnd_A-y_A+1];
        handles.XY_list_A(handles.smallcount_A,3)=timestamp_A;
    end
    
    % ANIMAL B FIX 2
    if x_B==1 && y_B==1;
        try
            handles.XY_list_B(handles.smallcount_B,1)= handles.XY_list_B(handles.smallcount_B-1,1);
            handles.XY_list_B(handles.smallcount_B,2)= handles.XY_list_B(handles.smallcount_B-1,2);
        catch
            handles.XY_list_B(handles.smallcount_B,1)= 1;
            handles.XY_list_B(handles.smallcount_B,2)= 1;
        end
        handles.XY_list_B(handles.smallcount_B,3)= timestamp_B;
    else
        handles.XY_list_B(handles.smallcount_B,1)=[x_B];
        handles.XY_list_B(handles.smallcount_B,2)=[handles.yEnd_B-y_B+1];
        handles.XY_list_B(handles.smallcount_B,3)=timestamp_B;
    end
    %}
    
    % ***_________ Want to keep this (these are the timing pulses) _________***
    % -------------------------------------------------------------------------
    %send the synchronization pulse
    
    if verbose == 1, disp('Sending synchronization pulse.'), end;
    pulseTimeOnset = unixtime(clock);
    StimTriggerAct('Pin35On');
    pause(2)
    StimTriggerAct('Pin35Off');
    
    %check for the proper directory
    today=date;
    %     if ~exist([handles.FILESAVEDIRECTORY_A filesep handles.ANIMAL_A_NUMBER '_' (today)],'dir');
    %         mkdir([handles.FILESAVEDIRECTORY_A filesep handles.ANIMAL_A_NUMBER '_' (today)])
    %     end
    %     if ~exist([handles.FILESAVEDIRECTORY_B filesep handles.ANIMAL_B_NUMBER '_' (today)],'dir');
    %         mkdir([handles.FILESAVEDIRECTORY_B filesep handles.ANIMAL_B_NUMBER '_' (today)])
    %     end
    
    %%%   don't need the XY list
    %     XY_list_A=[]; XY_list_B=[]; saveframe={};
    %     XY_list_A     =handles.XY_list_A;
    %     XY_list_B     =handles.XY_list_B;
    %%% Keep the saveframe
    try
        saveframe     = handles.saveframe;
    catch
        disp('broke')
        keyboard
    end
    %save the data
    %save([handles.FILESAVEDIRECTORY filesep handles.ANIMAL_NUMBER '_' (today) filesep handles.ANIMAL_NUMBER '_',handles.TODAYS_DATE,'_',num2str(handles.count)],...
    %   'XY_list','saveframe','pulseTimeOnset');
    % the save statement can omit XY list
    %    save ([handles.FILESAVEDIRECTORY_A filesep handles.ANIMAL_A_NUMBER '_' num2str(handles.count_A)],'XY_list_A','saveframe','pulseTimeOnset');
    %    save ([handles.FILESAVEDIRECTORY_B filesep handles.ANIMAL_B_NUMBER '_' num2str(handles.count_B)],'XY_list_B','saveframe','pulseTimeOnset');
    %     t01 = toc(t0);
    %     fprintf('%.3f seconds until now...\n',t01);
    %     t02 = tic;
    % if nAnimals == 1
    %     save ([handles.FILESAVEDIRECTORY filesep handles.ANIMAL_NUMBER ...
    %         '_' num2str(handles.count_A)],'saveframe','pulseTimeOnset');
    % elseif nAnimals == 2
    try
        save ([handles.FILESAVEDIRECTORY filesep handles.ANIMAL_NUMBER '_' num2str(handles.count_A)],...
            'saveframe','pulseTimeOnset','-v7.3');
    catch
        disp('here')
        keyboard;
    end
    % end
    %     save ([handles.FILESAVEDIRECTORY_B filesep handles.ANIMAL_B_NUMBER '_' num2str(handles.count_B)],'saveframe','pulseTimeOnset');
    %     t03 = toc(t02);
    
    flushdata(handles.vid,'triggers');
    
    %%%% RESET TRACKER - DO NEED THIS! I put this back
    % Not having this line makes each saved file contain all previous files
    % in one run -> crash!
    handles.resetTracker=1;
    
    %     t1 = toc(t0);
    if verbose == 1, fprintf('Exiting save..\n\n'), end;
else
    handles.imcount = handles.imcount + 1;
    % REGULAR TIMESTEP
    if verbose == 1, fprintf('Beginning regular timestep.\n'), end
    
    disp(['Frame snapped at: ' datestr(now,'dd-mmm-yyyy HH:MM:SS.FFF') '.']);
    %     clock
    
    %%% shouldn't need this:
    %{
% for some reason, the values here are not all round integers, some have
% decimals. this is a quick fix...
handles.yStart_A    = round(handles.yStart_A);
handles.yStart_B    = round(handles.yStart_B);
handles.yEnd_A      = round(handles.yEnd_A);
handles.yEnd_B      = round(handles.yEnd_B);
handles.xStart_A    = round(handles.xStart_A);
handles.xStart_B    = round(handles.xStart_B);
handles.xEnd_A      = round(handles.xEnd_A);
handles.xEnd_B      = round(handles.xEnd_B);
    %}
    
    
    % Only need the part where a frame is acquired:
    %{
% HERE CHANGE SO IT DEFINES 2 AREAS TO LOOK INTO
    all_yStart  = min(handles.yStart_A,handles.yStart_B);
    all_yEnd    = max(handles.yEnd_A,handles.yEnd_B);
    all_xStart  = min(handles.xStart_A,handles.xStart_B);
    all_xEnd    = max(handles.xEnd_A,handles.xEnd_B);
    %}
    trigger (handles.vid)
    %     pause(0.2)
    datT            = clock;
    
    verb = 0;
    while get(handles.vid,'FramesAvailable')<1
        if verb==1
            tic
            disp('Frames unavailable. Waiting before retry.')
            toc
        end
    end
    
    svdata          = getdata(handles.vid,1);
    
    % this added to check camera field of view
    if handles.viewIm
        if handles.imcount > 10 % every 10 frames
            data = svdata(:,:,3);
            img = imresize(data,0.2,'bilinear');
            figure(800);
            imshow(img);
            handles.imcount  = 0;
            pause(.1);
        end
    end
    
    % I think these variable are only needed for the XY thing, so can comment
    % out:
    %{
    data            = (svdata(:,:,3));
    %data            = rgb2gray(data);
    data_all        = data(all_yStart:all_yEnd,all_xStart:all_xEnd);
    data_A          = data(handles.yStart_A:handles.yEnd_A,handles.xStart_A:handles.xEnd_A);
    data_B          = data(handles.yStart_B:handles.yEnd_B,handles.xStart_B:handles.xEnd_B);
    data_all        = imresize(data_all,0.2,'bilinear');
    data_A          = imresize(data_A,0.2,'bilinear');
    data_B          = imresize(data_B,0.2,'bilinear');
    %}
    
    %%% Can remove this entirely:
    %{
    % if filters not defined, make them
    if ~isfield(handles,'filters_A') || ~isfield(handles,'filters_B');
        figure(78);
        imshow(data_A);
        txtA=text(0, 150, 'Click the center of each bulb for animal A (LEFT)',...
            'Fontsize',18,'Color','green');
        [bulbsX_A bulbsY_A] = ginput(2);
        delete(txtA);
        close(figure(78));
        figure(79);
        imshow(data_B);
        txtB=text(0, 150, 'Click the center of each bulb for animal B (RIGHT)',...
            'Fontsize',18,'Color','green');
        [bulbsX_B bulbsY_B] = ginput(2);
        delete(txtB);
        close (figure (79));
        SEP_A = pdist([bulbsX_A bulbsY_A]);
        SEP_B = pdist([bulbsX_B bulbsY_B]);
        RAD = 4;
        DA  = 18;
        [FILTERS_A,X_A,Y_A] = dual_led_filters(DA, SEP_A, RAD);
        [FILTERS_B,X_B,Y_B] = dual_led_filters(DA, SEP_B, RAD);
        handles.filters_A = FILTERS_A;
        handles.filters_B = FILTERS_B;
    end

    % find x and y coordinates of IR bulbs for each animal
    [x_A,y_A,timestamp_A] = alpha_v3 (data_A,handles.filters_A);
    [x_B,y_B,timestamp_B] = alpha_v3 (data_B,handles.filters_B);
 
    
    % ANIMAL A FIX 1
    if length(x_A)~=1 || length (y_A)~=1;
        disp('XY problem detected and trying to fix');
        
        if length(x_A)>1;
           x_A = x_A(1);
        end
        if length(y_A)>1;
            y_A = y_A(1);
        end
        if isempty(x_A);
            x_A = 1;
        end
        if isempty(y_A);
            y_A = 1;
        end
    end
    
    % ANIMAL B FIX 1
    if length(x_B)~=1 || length (y_B)~=1;
        disp('XY problem detected and trying to fix');
        
        if length(x_B)>1;
            x_B = x_B(1);
        end
        if length(y_B)>1;
            y_B = y_B(1);
        end
        if isempty(x_B);
            x_B = 1;
        end
        if isempty(y_B);
            y_B = 1;
        end
    end
    
    % ANIMAL A FIX 2 and SAVE DATA
    if x_A==1 && y_A==1;
        
        if handles.count_A == 1;
            error('VIDACQUISITION:timeVid_v2:nogo','Can not find a light on the first frame. Check the bulb and restart.');
        end
        
        sc = handles.smallcount_A-1;
        if sc == 0;
            % this is an annoying scenario in which the first frame of a
            % batch finds XY as [1,1]. Thus, we can't take the last known
            % position because it doesn't exist. So, we'll try grabbing
            % another frame, and using those coordinates, even if they
            % still return [1,1]. These erronious points can be dealt with
            % after the recording has been completed.
            [x_A,y_A,timestamp_A] = alpha_v3 (data_A,handles.filters_A);
            handles.XY_list_A(handles.smallcount_A,1)   = [x_A];
            handles.XY_list_A(handles.smallcount_A,2)   = [y_A];%[handles.yEnd-y+1];
            handles.XY_list_A(handles.smallcount_A,3)   = timestamp_A;
            
        else
            % this is the better scenario for dealing with coordinates of
            % [1,1].
            handles.XY_list_A(handles.smallcount_A,1)   = handles.XY_list_A(sc,1);
            handles.XY_list_A(handles.smallcount_A,2)   = handles.XY_list_A(sc,2);
            handles.XY_list_A(handles.smallcount_A,3)   = timestamp_A;
            
        end

    else
        handles.XY_list_A(handles.smallcount_A,1)   = x_A;
        handles.XY_list_A(handles.smallcount_A,2)   = y_A;%[handles.yEnd-y+1];
        handles.XY_list_A(handles.smallcount_A,3)   = timestamp_A;
    end
    
    % ANIMAL B FIX 2 and SAVE DATA
    if x_B==1 && y_B==1;
        
        if handles.count_B == 1;
            error('VIDACQUISITION:timeVid_v2:nogo','Can not find a light on the first frame. Check the bulb and restart.');
        end
        
        sc = handles.smallcount_B-1;
        if sc == 0;
            % this is an annoying scenario in which the first frame of a
            % batch finds XY as [1,1]. Thus, we can't take the last known
            % position because it doesn't exist. So, we'll try grabbing
            % another frame, and using those coordinates, even if they
            % still return [1,1]. These erronious points can be dealt with
            % after the recording has been completed.
            [x_B,y_B,timestamp_B] = alpha_v3 (data_B,handles.filters_B);
            handles.XY_list_B(handles.smallcount_B,1)   = [x_B];
            handles.XY_list_B(handles.smallcount_B,2)   = [y_B];%[handles.yEnd-y+1];
            handles.XY_list_B(handles.smallcount_B,3)   = timestamp_B;
            
        else
            % this is the better scenario for dealing with coordinates of
            % [1,1].
            handles.XY_list_B(handles.smallcount_B,1)   = handles.XY_list_B(sc,1);
            handles.XY_list_B(handles.smallcount_B,2)   = handles.XY_list_B(sc,2);
            handles.XY_list_B(handles.smallcount_B,3)   = timestamp_B;
            
        end

    else
        handles.XY_list_B(handles.smallcount_B,1)   = x_B;
        handles.XY_list_B(handles.smallcount_B,2)   = y_B;%[handles.yEnd-y+1];
        handles.XY_list_B(handles.smallcount_B,3)   = timestamp_B;
    end

end
    %}
    
    % This is the part I need to change to get a higher frame rate.
    % Right now the save happens every time count_master is a multiple of 10.
    % This is an imprecise way of timing the save process - resulting in
    % irregular time intervals between frames. I have to find a better way to
    % do this, based on timing rather than a counter.
    % I'll try to turn up the frame rate to be 1 frame/sec initially. If this
    % is feasible I can try to increase to somewhere between 1 and 30
    % frames/sec, depending on how well the timing works.
    %
    % 11/23/15 3:00 PM FIX FOUND:
    %  the time interval is now precisely what indicated by the 'Period'
    %  property of the timer object defined in videotimer_v3_2animals.m
    %
    % There is still a larger delay (~6 sec) when a save event occurs (usually
    % we set this to every 1/2 hour).
    %
    %
    %%%% **********************************************************************
    %%%% SAVEFRAME IS HERE
    %%%% **********************************************************************
    %save a black and white version of the video frame every 5-6 seconds or so
    
    if verbose == 1, disp('storing image data'), end
    
    handles.framecount                              = handles.framecount+1;
    tempframe                                       = [];
    tempframe                                       = rgb2gray(svdata);
    handles.saveframe.image{handles.framecount}     = imresize(tempframe, 0.2);
    save_time                                       = [];
    
    if verbose == 1;
        disp(['Frame ' num2str(handles.framecount) ' timestamp is: ' ...
            datestr(datT, 'dd-mmm-yyyy HH:MM:SS.FFF')]);
    end
    
    save_time                                       = unixtime(datT);
    handles.saveframe.time(handles.framecount)      = save_time;
    
    uploadchecker_pi(save_time,handles.checkerDir);
    
    
end

%%% Can get rid of this:
%{
%--------------------------------------------------------------------------
%%% This is for the self-updating rat tracker 3000 plot
%}
%%%% PUTTING THIS BACK up there
if handles.resetTracker     == 1;
    
    handles.saveframe.image    = {};
    handles.saveframe.time     = [];
    
    handles.framecount         = 0;
    handles.resetTracker       = 0;
end
%{
    try
    delete(handles.ratTracker3000_A);
    catch
        keyboard
    end
    try
    delete(handles.ratTracker3000_B);
    catch
        keyboard
    end
    handles.resetTracker       = 0;

else
    
    if rem(handles.count_master,15) == 0;
        
        mins_A=handles.smallcount_A-500;
        if mins_A<1;
            mins_A=1;
        end
        mins_B=handles.smallcount_B-500;
        if mins_B<1;
            mins_B=1;
        end
        
        if verbose == 1, disp('about to plot.'), end;
        % Animal A (left) RatTracker
        handles.ratTracker3000_A = figure(3);
        %data2=svdata(handles.yStart:handles.yEnd,handles.xStart:handles.xEnd,3);
        imshow(data_A), hold on;
        plot(handles.XY_list_A(mins_A:handles.smallcount_A,1),handles.XY_list_A(mins_A:handles.smallcount_A,2)); hold on
        try
            delete (handles.f_A);
        end
        handles.f_A = scatter(x_A,y_A,'rx');
        set(gcf, 'MenuBar', 'none','Position',[200 500 400 400]);
        set(gcf, 'ToolBar', 'none');
        set(gca,'XTick',[],'YTick',[]);
        title('Rat Tracker 3000 - animal A (left)');
        
        % Animal B (right) RatTracker
        handles.ratTracker3000_B = figure(4);
        imshow(data_B), hold on;
        plot(handles.XY_list_B(mins_B:handles.smallcount_B,1),handles.XY_list_B(mins_B:handles.smallcount_B,2)); hold on
        try
            delete (handles.f_B);
        end
        handles.f_B = scatter(x_B,y_B,'rx');
        set(gcf, 'MenuBar', 'none','Position',[700 500 400 400]);
        set(gcf, 'ToolBar', 'none');
        set(gca,'XTick',[],'YTick',[]);
        title('Rat Tracker 3000 - animal B (right)');
    end
    
end
if verbose ==1,disp('done plotting.'), end;
%}
