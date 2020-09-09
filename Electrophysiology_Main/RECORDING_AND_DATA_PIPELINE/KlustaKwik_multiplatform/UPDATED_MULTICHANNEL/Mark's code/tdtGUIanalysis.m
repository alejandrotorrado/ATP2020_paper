function tdtGUIanalysis()
% tdtGUIanalysis - graphical user interface to multichannel analysis, to be reincorporated into single tdtGUI application
% TODO (tdtGUIanalysis) fix infinite loop race condition -> no repeated variables (find ' = ' and 'function'), especially loop indices
% TODO (tdtGUIanalysis) do all graph updates in cbTimer -> set flags elsewhere, eventually allow cancellations (activity list with buttons)

  % initialize data
  iNerves    = 2;
  iProbes    = 4;
  iChannels  = iProbes*16;
  iChannel   = [];
  dRawData   = cell(1,iChannels);
  zRawData   = 0;
  dDecData   = [];
  iDecData   = 0;
  zDecData   = 0;
  dSpikes    = [];
  tSpikes    = [];
  sSpikes    = [];
  hSpikes    = [];
  ySpikes    = ones(1,iChannels)*400;
  dNerves    = cell(1,iNerves);
  zNerves    = 0;
  dDecNerves = [];
  iDecNerves = 0;
  dAllNerves = [];
  dCycles    = [];
  tCycles    = [];
  yCycles    = ones(1,iNerves)*2000;
  pDown      = [];
  tStart     = clock;
  fActive    = ['Data',filesep,datestr(date,29)];
  if exist(fActive,'file')
    i = 1;
    while exist([fActive,'-',int2str(i)],'file')
      i = i+1;
    end
    fActive = [fActive,'-',int2str(i)];
  end
  mkdir(fActive);
  fActive    = [pwd,filesep,fActive,filesep];
  tdtZBus    = [];
  tdtRX5     = [];
  tdtRP2     = 1;

  % initialize GUI
  setGUIDefaults();
  hWindow                  = figure( ...
    'Name',                  'Multichannel Data Acquisition and Analysis', ...
    'Color',                 'black', ...
    'CloseRequestFcn',       @cbClose, ...
    'KeyPressFcn',           @cbKey, ...
    'WindowButtonDownFcn',   @cbMouseDown, ...
    'WindowButtonMotionFcn', @cbMouseMove, ...
    'WindowButtonUpFcn',     @cbMouseUp);
  hTimer                   = timer( ...
    'Period',                0.1, ...
    'UserData',              0.25, ...
    'ExecutionMode',         'fixedDelay', ...
    'TimerFcn',              @cbTimer);

  % initialize toolbar
  hToolbar                 = uitoolbar();
  hToolFolderOpen          = uipushtool( ...
    'CData',                 getIcon('folder-open'), ...
    'TooltipString',         'Choose Folder to Open', ...
    'ClickedCallback',       @cbSetFolder);
  hToolFolderSave          = uipushtool( ...
    'CData',                 getIcon('folder-save'), ...
    'TooltipString',         'Choose Folder for Saving', ...
    'ClickedCallback',       @cbSetFolder);
  % TODO (hToolPrint,hToolPrintPreview)
  hToolRecord              = uitoggletool( ...
    'Separator',             'on', ...
    'CData',                 getIcon('status-record'), ...
    'TooltipString',         'Record', ...
    'ClickedCallback',       @cbSetState);
  hToolPlay                = uitoggletool( ...
    'CData',                 getIcon('status-play'), ...
    'TooltipString',         'Play', ...
    'ClickedCallback',       @cbSetState);
  hToolIdle                = uitoggletool( ...
    'State',                 'on', ...
    'CData',                 getIcon('status-idle'), ...
    'TooltipString',         'Idle', ...
    'ClickedCallback',       @cbSetState);
  hToolProbe1              = uitoggletool( ...
    'Separator',             'on', ...
    'CData',                 getIcon('probe-1'), ...
    'TooltipString',         '1 Multichannel Probe', ...
    'ClickedCallback',       @cbSetProbeNumber);
  hToolProbe2              = uitoggletool( ...
    'CData',                 getIcon('probe-2'), ...
    'TooltipString',         '2 Multichannel Probes', ...
    'ClickedCallback',       @cbSetProbeNumber);
  hToolProbe3              = uitoggletool( ...
    'CData',                 getIcon('probe-3'), ...
    'TooltipString',         '3 Multichannel Probes', ...
    'ClickedCallback',       @cbSetProbeNumber);
  hToolProbe4              = uitoggletool( ...
    'State',                 'on', ...
    'CData',                 getIcon('probe-4'), ...
    'TooltipString',         '4 Multichannel Probes', ...
    'ClickedCallback',       @cbSetProbeNumber);
  hToolNerve1              = uitoggletool( ...
    'Separator',             'on', ...
    'CData',                 getIcon('nerve-1'), ...
    'TooltipString',         '1 Nerve Electrode', ...
    'ClickedCallback',       @cbSetNerveNumber);
  hToolNerve2              = uitoggletool( ...
    'State',                 'on', ...
    'CData',                 getIcon('nerve-2'), ...
    'TooltipString',         '2 Nerve Electrodes', ...
    'ClickedCallback',       @cbSetNerveNumber);
  hToolAnalysis            = uitoggletool( ...
    'Separator',             'on', ...
    'State',                 'on', ...
    'CData',                 getIcon('tab-analysis'), ...
    'TooltipString',         'Analysis', ...
    'ClickedCallback',       @cbSetPositions);
  hToolStatus              = uitoggletool( ...
    'State',                 'on', ...
    'CData',                 getIcon('tab-status'), ...
    'TooltipString',         'Status', ...
    'ClickedCallback',       @cbSetPositions);

  % initialize viewer panel
  hBorderSettings          = axes();
  hTextClock               = uicontrol( ...
    'Style',                 'text', ...
    'String',                datestr(clock), ...
    'HorizontalAlignment',   'center');
  hTextState               = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'right');
  hTextStateTime           = uicontrol( ...
    'Style',                 'text');
  hLabelFile               = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Current File:', ...
    'Visible',               'off', ...
    'HorizontalAlignment',   'right');
  hTextFile                = uicontrol( ...
    'Style',                 'text', ...
    'UserData',              0, ...
    'HorizontalAlignment',   'right');
  hLabelGraphSettings      = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Graph Settings', ...
    'HorizontalAlignment',   'center');
  hLabelTimeScale          = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Time:', ...
    'HorizontalAlignment',   'right');
  hMenuTimeScale           = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'1 Second';'2 Seconds';'5 Seconds';'10 Seconds';'15 Seconds';'20 Seconds';'30 Seconds';'1 Minute';'2 Minutes';'5 Minutes';'10 Minutes';'15 Minutes';'20 Minutes';'30 Minutes';'1 Hour';'2 Hours';'All Data'}, ...
    'UserData',              [1,2,5,10,15,20,30,60,120,300,600,900,1200,1800,3600,7200,0], ...
    'Value',                 10, ...
    'Callback',              @cbSetTimeScales);
  hLabelNerveVoltScale     = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Nerve:', ...
    'HorizontalAlignment',   'right');
  hMenuNerveVoltScale      = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'200 mV';'400 mV';'600 mV';'800 mV';'1 V';'1.5 V';'2 V';'2.5 V';'3 V';'4 V';'5 V';'10 V'}, ...
    'UserData',              [2000,4000,6000,8000,10000,15000,20000,25000,30000,40000,50000,100000], ...
    'Value',                 5, ...
    'Callback',              @cbSetNerveVoltScales);
  hLabelDataVoltScale      = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Probe:', ...
    'HorizontalAlignment',   'right');
  hMenuDataVoltScale       = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'25 uV';'50 uV';'75 uV';'100 uV';'150 uV';'200 uV';'300 uV';'400 uV';'500 uV';'750 uV';'1 mV';'3.3 mV'}, ...
    'UserData',              [250,500,750,1000,1500,2000,3000,4000,5000,7500,10000,33000], ...
    'Value',                 4, ...
    'Callback',              @cbSetDataVoltScales);
  hLabelView               = uicontrol( ...
    'Style',                 'text', ...
    'String',                'View:', ...
    'HorizontalAlignment',   'right');
  hMenuView                = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'All';'1-16';'17-32';'33-48';'49-64'}, ...
    'Callback',              @cbSetPositions);
  hSliderTime              = uicontrol( ...
    'Style',                 'slider', ...
    'BackgroundColor',       [.7,.7,.7], ...
    'Callback',              @cbSetTime);
  hMenuTime                = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'1/8 x';'1/4 x';'1/2 x';'Real Time';'2 x';'4 x';'8 x'}, ...
    'Value',                 4);
  hLabelNerve              = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Nerve(s)');
  % TODO (hGraphNerve) label bursts (highlight or underline) -> interface to ignore/split
  hGraphNerve              = axes();
  hGraphNerveLine          = ones(1,iNerves);
  for i = 1:iNerves
    hGraphNerveLine(i)     = line( ...
      'Color',               getColor(i));
  end
  hLabelData               = ones(1,iChannels);
  hGraphData               = ones(1,iChannels);
  hGraphDataLine           = ones(1,iChannels);
  hGraphSpikes             = ones(1,iChannels);
  hGraphSpikesLine         = ones(1,iChannels);
  for i = 1:iChannels
    hLabelData(i)          = uicontrol( ...
      'Style',               'text', ...
      'String',              i);
    hGraphData(i)          = axes( ...
      'UserData',            i, ...
      'ButtonDownFcn',       @cbSetChannel);
    hGraphDataLine(i)      = line( ...
      'Color',               getColor(i));
    hGraphSpikes(i)        = axes( ...
      'UserData',            i, ...
      'XLim',                [1,62], ...
      'ButtonDownFcn',       @cbSetChannel);
    hGraphSpikesLine(i)    = line( ...
      'Color',               getColor(i));
  end
  hGraphScope              = axes( ...
    'Color',                 'none', ...
    'HitTest',               'off');
  hGraphScopeLine          = line( ...
    'XData',                 [0,0], ...
    'YData',                 [0,1]);

  % initialize analysis panel
  hLabelCycle              = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Cycles');
  hTextCycle               = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'right');
  % TODO (hMenuCycle) raise/lower with mouse and save/edit with button (disable while recording), change auto to reset to default
  hMenuCycle               = uicontextmenu();
  hMenuCycleSave           = uimenu(hMenuCycle, ...
    'Label',                 'Save Cycles', ...
    'Callback',              @cbSetThreshold);
  hMenuCycleAuto           = uimenu(hMenuCycle, ...
    'Label',                 'Auto Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuCycleRaise          = uimenu(hMenuCycle, ...
    'Label',                 'Raise Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuCycleLower          = uimenu(hMenuCycle, ...
    'Label',                 'Lower Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuCycle2              = uicontextmenu();
  hMenuCycleEdit           = uimenu(hMenuCycle2, ...
    'Label',                 'Edit Cycle Threshold', ...
    'Callback',              @cbSetThreshold);
  % TODO (hGraphCycle) stacked cycles OR aligned with threshold graph
  hGraphCycle              = axes( ...
    'UIContextMenu',         hMenuCycle);
  hGraphCycleLine          = [];
  hLabelThresh             = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Spike Threshold');
  hTextThresh              = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'right');
  % TODO (hMenuThresh) raise/lower with mouse and save/edit (disable while recording),auto,all,ignore with buttons
  hMenuThresh              = uicontextmenu();
  hMenuThreshSave          = uimenu(hMenuThresh, ...
    'Label',                 'Save Spikes', ...
    'Callback',              @cbSetThreshold);
  hMenuThreshAutoAll       = uimenu(hMenuThresh, ...
    'Label',                 'Auto Threshold All', ...
    'Callback',              @cbSetThreshold);
  hMenuThreshAuto          = uimenu(hMenuThresh, ...
    'Label',                 'Auto Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuThreshRaise         = uimenu(hMenuThresh, ...
    'Label',                 'Raise Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuThreshLower         = uimenu(hMenuThresh, ...
    'Label',                 'Lower Threshold', ...
    'Callback',              @cbSetThreshold);
  hMenuThreshIgnore        = uimenu(hMenuThresh, ...
    'Label',                 'Ignore Channel', ...
    'Callback',              @cbSetThreshold);
  hMenuThresh2             = uicontextmenu();
  hMenuThreshEdit          = uimenu(hMenuThresh2, ...
    'Label',                 'Edit Spike Thresholds', ...
    'Callback',              @cbSetThreshold);
  hGraphThresh             = axes( ...
    'UIContextMenu',         hMenuThresh);
  hGraphThreshLine         = line();
  hGraphThreshLimitLine    = line( ...
    'Color',                 [.5,.5,.5]);
  hGraphThreshLimitLineNeg = line( ...
    'Color',                 [.5,.5,.5]);
  hLabelSort               = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Spike Sort');
  hGraphSort               = axes( ...
    'UserData',              1, ...
    'XLim',                  [1,62]);
  hGraphSortLine           = ones(1,30);
  for i = 1:30
    hGraphSortLine(i) = line();
  end
  hCheckSortEasy           = uicontrol( ...
    'Style',                 'checkbox', ...
    'String',                'Filter Length and Amplitude', ...
    'Callback',              @cbSetSort);
  hCheckSortLow            = uicontrol( ...
    'Style',                 'checkbox', ...
    'String',                'Filter Low Amplitude Harshly', ...
    'Value',                 0, ...
    'Callback',              @cbSetSort);
  hCheckSortOutliers       = uicontrol( ...
    'Style',                 'checkbox', ...
    'String',                'Discard Outliers', ...
    'Callback',              @cbSetSort);
  hMenuSortDirection       = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'Directions combined','Rising only','Falling only','Directions split'}, ...
    'Callback',              @cbSetSort);
  hGraphAux                = axes();
  hTextAux                 = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'right');
  hMenuAux                 = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'Cycle Triggered Histogram';'Spike Triggered Histogram';'Interspike Intervals';'Interspike Intervals Histogram';'Firing Rate Histogram'}, ...
    'Callback',              @cbSetAux);
  hMenuAuxCTH              = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'Cycles from Nerve 1','Cycles from Nerve 2'}, ...
    'Callback',              @cbSetAux);
  hMenuOutput              = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'Output Spike Frequencies','Output Spike Times'});
  hButtonOutput            = uicontrol( ...
    'String',                'Print Output', ...
    'ForegroundColor',       'white', ...
    'Callback',              @cbSetOutput);
  setChannel(1);
  cbSetNerveVoltScales(hMenuNerveVoltScale,[]);
  cbSetDataVoltScales(hMenuDataVoltScale,[]);

  % initialize status panel
  hBorderRX5               = axes();
  hLabelRX5                = uicontrol( ...
    'Style',                 'text', ...
    'String',                'RX5 Status', ...
    'HorizontalAlignment',   'center');
  hTextRX5Status           = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'center');
  hTextRX5Data             = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'center');
  hLabelRX5Filter          = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Band Pass Filter', ...
    'HorizontalAlignment',   'center');
  hMenuRX5FilterHighPass   = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'.1 Hz';'1 Hz';'10 Hz';'100 Hz';'200 Hz';'300 Hz';'400 Hz'}, ...
    'UserData',              [.1,1,10,100,200,300,400], ...
    'Value',                 6, ...
    'Callback',              @cbSetFilter);
  hLabelRX5FilterDash      = uicontrol( ...
    'Style',                 'text', ...
    'String',                '-', ...
    'HorizontalAlignment',   'center');
  hMenuRX5FilterLowPass    = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                {'500 Hz';'1 kHz';'2 kHz';'3 kHz';'4 kHz';'5 kHz';'10 kHz'}, ...
    'UserData',              [500,1000,2000,3000,4000,5000,10000], ...
    'Value',                 6, ...
    'Callback',              @cbSetFilter);
  hBorderRP2               = axes();
  hLabelRP2                = uicontrol( ...
    'Style',                 'text', ...
    'String',                'RP2 Status', ...
    'HorizontalAlignment',   'center');
  hTextRP2Status           = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'center');
  hTextRP2Data             = uicontrol( ...
    'Style',                 'text', ...
    'HorizontalAlignment',   'center');
  hBorderUsages            = axes();
  hLabelUsages             = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Processor Usages', ...
    'HorizontalAlignment',   'center');
  hGraphUsages             = axes( ...
    'CLimMode',              'auto', ...
    'XLim',                  [.5,3.5], ...
    'XTick',                 [1,2,3], ...
    'XTickLabel',            {'RX5-Main','RX5-Aux','RP2'}, ...
    'YLim',                  [0,100]);
  hGraphUsageRX5M          = patch([0.7,1.3,1.3,0.7],[0,0,0,0],[1,1,1]);
  hGraphUsageRX5A          = patch([1.7,2.3,2.3,1.7],[0,0,0,0],[1,1,1]);
  hGraphUsageRP2           = patch([2.7,3.3,3.3,2.7],[0,0,0,0],[1,1,1]);
  set(get(hGraphUsages,'Children'), ...
    'EdgeColor',             'none', ...
    'FaceColor',             'interp', ...
    'FaceVertexCData',       [0,0,0;0,0,0;1,1,1;1,1,1]);
  hBorderEvents            = axes();
  hLabelEvents             = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Events', ...
    'HorizontalAlignment',   'center');
  hButtonEvent             = uicontrol( ...
    'String',                'Add', ...
    'ForegroundColor',       'black', ...
    'Callback',              @cbSetEvent);
  hMenuEvents              = uicontrol( ...
    'Style',                 'popupmenu', ...
    'String',                'Other');
  hTextEvents              = uicontrol( ...
    'Style',                 'edit');
  hBorderMessages          = axes();
  hLabelMessages           = uicontrol( ...
    'Style',                 'text', ...
    'String',                'Messages', ...
    'HorizontalAlignment',   'center');
  hTextMessages            = uicontrol( ...
    'Style',                 'edit', ...
    'Max',                   2, ...
    'Enable',                'inactive');

  % make the window visible and start timer
  cbSetPositions([],[]);
  setMessage('Program started');
  setMessage(['Using folder: ',fActive]);
  setTDTDefaults(1);
  set(hWindow,'Visible','on');
  setCameras(hWindow);
  start(hTimer);

  % callback for close request
  function cbClose(source,event)
    bClose = false;
    if ~strcmp(getState(),'recording')
      bClose = true;
    elseif strcmp(questdlg('Stop recording and exit?','Close ...','OK','Cancel','OK'),'OK')
      set(hToolIdle,'State','on');
      cbSetState(hToolIdle,event);
      bClose = true;
    end
    if bClose
      set(hWindow,'Visible','off');
      stop(hTimer);
      setMessage('Program closed');
      if ~exist([fActive,'messages.mat'],'file')
        try
          rmdir(fActive);
        catch
          lasterror
        end
      end
      delete(hTimer);
      delete(source);
      clear all;
    end
  end

  % callback for key press
  % TODO (cbKey) additional keys for navigating timeline, especially during replay
  function cbKey(source,event)
    k = lower(event.Key);
    m = event.Modifier;
    if (strcmp(k,'w') && length(m) == 1 && strcmp(m{1},'control')) || strcmp(k,'escape')
      close(gcf);
    elseif length(k) == 1 && length(m) == 1 && strcmp(m{1},'control')
      if strcmp(k,'o') && strcmp(get(hToolFolderOpen,'Enable'),'on')
        cbSetFolder(hToolFolderOpen,event);
      elseif strcmp(k,'s') && strcmp(get(hToolFolderSave,'Enable'),'on')
        cbSetFolder(hToolFolderSave,event);
      elseif strcmp(k,'r') && strcmp(get(hToolRecord,'Enable'),'on') && strcmp(get(hToolRecord,'State'),'off')
        set(hToolRecord,'State','on');
        cbSetState(hToolRecord,event);
      elseif strcmp(k,'p') && strcmp(get(hToolPlay,'Enable'),'on') && strcmp(get(hToolPlay,'State'),'off')
        set(hToolPlay,'State','on');
        cbSetState(hToolPlay,event);
      elseif strcmp(k,'i') && strcmp(get(hToolIdle,'Enable'),'on') && strcmp(get(hToolIdle,'State'),'off')
        set(hToolIdle,'State','on');
        cbSetState(hToolIdle,event);
      end
    else
      v = iChannel;
      if strcmp(k,'pageup') || strcmp(event.Character,'+')
        v = min(v+4,iChannels);
      elseif strcmp(k,'uparrow')
        v = min(v+1,iChannels);
      elseif strcmp(k,'downarrow')
        v = max(v-1,1);
      elseif strcmp(k,'pagedown') || strcmp(event.Character,'-')
        v = max(v-4,1);
      end
      if v ~= iChannel
        setChannel(v);
      end
    end
  end

  % callback for mouse click
  % TODO (cbMouseDown) if in hGraphCycle, edit cycle threshold
  % TODO (cbMouseDown) if in hGraphThresh, edit spike threshold
  % TODO (cbMouseDown) if in hGraphSort, edit time/voltage windows
  % TODO (cbMouseDown) if in hGraphAux, edit pre/E/I/post lines (initialize to 10% height, 80% time, 10% height)
  function cbMouseDown(source,event)
    pDown = get(source,'CurrentPoint');
  end

  % callback for mouse move or drag
  % TODO (cbMouseMove) if in hGraphNerve, highlight nearby cycle in hGraphCycle
  % TODO (cbMouseMove) if in hGraphThresh, highlight nearby spike in hGraphSort
  function cbMouseMove(source,event)
    if isempty(pDown)
    else
    end
  end

  % callback for mouse release
  % TODO (cbMouseUp) update with results of edit
  function cbMouseUp(source,event)
    pDown = [];
  end

  % callback for auxiliary graph menu
  function cbSetAux(source,event)
    if get(hMenuAux,'Value') == 1 && iNerves > 1
      set(hMenuAuxCTH,'Visible','on');
    else
      set(hMenuAuxCTH,'Visible','off');
    end
    t = get(hSliderTime,'Value');
    setAux(getTimeMinimum(t),t);
  end

  % callback for clicking on a channel
  function cbSetChannel(source,event)
    setChannel(get(source,'UserData'));
  end

  % callback for data volt scale menu
  function cbSetDataVoltScales(source,event)
    n = get(source,'UserData');
    n = n(get(source,'Value'));
    set([hGraphData,hGraphSpikes,hGraphThresh,hGraphSort],'YLim',[-n,n]);
    setCameras([hGraphData,hGraphSpikes,hGraphThresh,hGraphSort]);
  end

  % callback for adding event to message list
  function cbSetEvent(source,event)
    s = get(hMenuEvents,'String');
    if iscell(s)
      s = s{get(hMenuEvents,'Value')};
    end
    n = get(hTextEvents,'String');
    if ~isempty(n)
      s = [s,' - ',n];
    end
    setMessage(s);
    % TODO (cbSetEvent) add/edit list when replaying
    if strcmp(getState(),'recording')
      e = get(hTextEvents,'UserData');
      e(end+1).time = get(hSliderTime,'Value');
      e(end).string = s;
      set(hTextEvents,'UserData',e);
    end
  end

  % callback for multichannel band pass filter
  function cbSetFilter(source,event)
    if tdtIsReady()
      tdtSetFilterFreqs();
    end
  end

  % callback for folder tools
  function cbSetFolder(source,event)
    f = uigetdir(fActive);
    if f ~= 0
      if source == hToolFolderOpen
        if exist([f,filesep,'messages.mat'],'file')
          if ~exist([fActive,'messages.mat'],'file')
            try
              rmdir(fActive);
            catch
              lasterror
            end
          end
          fActive = [f,filesep];
          set([hToolRecord,hToolProbe1,hToolProbe2,hToolProbe3,hToolProbe4,hToolNerve1,hToolNerve2],'Enable','off');
          set(hTextRX5Status,'String','RX5 is offline.');
          set(hTextRP2Status,'String','RP2 is offline.');
          set([hGraphUsageRX5M,hGraphUsageRX5A,hGraphUsageRP2],'YData',[0,0,0,0]);
          set(hToolIdle,'State','on');
          cbSetState(hToolIdle,event);
          % TODO (cbSetFolder) load correctly if data.mat exists or multiple data-xxx.mat files
          m = [];
          load([fActive,'messages']);
          s = textscan(m,'%[^\n]');
          m = [];
          for s = s{:,1}'
            if (isempty(strfind(s{1},'Saved')) || strfind(s{1},'Saved') ~= 10) && (isempty(strfind(s{1},'Filled')) || strfind(s{1},'Filled') ~= 10)
              if isempty(m)
                m = s{1};
              else
                m = sprintf('%s\n%s',m,s{1});
              end
            end
          end
          set(hTextMessages,'String',m);
          tdtData = [];
          rawdata = [];
          nerve = [];
          d = dir([fActive,'data*.mat']);
          if ~isempty(d)
            load([fActive,d(1).name]);
          end
          tdtZBus = [];
          htn = hToolNerve1;
          for data = tdtData
            switch data{1}.name
              case 'Wave'
                tdtRX5 = data{1};
                set(hTextRX5Data,'String',sprintf('%d channels @ %d Hz',tdtRX5.channels,round(tdtRX5.samplingRate)));
                switch tdtRX5.channels/16
                  case 1
                    h = hToolProbe1;
                  case 2
                    h = hToolProbe2;
                  case 3
                    h = hToolProbe3;
                  case 4
                    h = hToolProbe4;
                end
                if strcmp(get(h,'State'),'off')
                  set(h,'State','on');
                  cbSetProbeNumber(h,event);
                end
                r = dir([fActive,'rawdata*.mat']);
                if ~isempty(r)
                  set(hSliderTime,'Min',0,'Max',length(r)*tdtRX5.bufferSize/tdtRX5.samplingRate,'Value',0);
                end
                decimated = [];
                if exist([f,filesep,'decimated.mat'],'file')
                  load([fActive,'decimated']);
                  iDecData = decimated.dataDecimation;
                  iDecNerves = decimated.nerveDecimation;
                  dDecData = decimated.data;
                  dDecNerves = decimated.nerves;
                  zDecData = length(dDecData{1});
                  set(hSliderTime,'Min',0,'Max',zDecData*.5*tdtRX5.blockSize/tdtRX5.samplingRate,'Value',0);
                else
                  dDecData = [];
                  iDecData = 0;
                  zDecData = 0;
                  dDecNerves = [];
                  iDecNerves = 0;
                end
              case 'Nrv1'
                tdtRP2 = data{1};
                set(hTextRP2Data,'String',sprintf('1 channel @ %d Hz',round(tdtRP2.samplingRate)));
              case 'Nrv2'
                htn = hToolNerve2;
                set(hTextRP2Data,'String',sprintf('2 channels @ %d Hz',round(tdtRP2.samplingRate)));
            end
          end
          if strcmp(get(htn,'State'),'off')
            set(htn,'State','on');
            cbSetNerveNumber(htn,event);
          end
          spikes = [];
          times = [];
          sortcodes = cell(1,iChannels);
          hSpikes = cell(1,iChannels);
          thresholds = ones(1,iChannels)*400;
          if exist([f,filesep,'spikes.mat'],'file')
            load([fActive,'spikes']);
            dSpikes = spikes;
            tSpikes = times;
            sSpikes = sortcodes;
            ySpikes = thresholds;
            set(hGraphThresh,'UIContextMenu',hMenuThresh2);
          else
            dSpikes = [];
            tSpikes = [];
            sSpikes = [];
            ySpikes = ones(1,iChannels)*400;
            set(hGraphThresh,'UIContextMenu',hMenuThresh);
          end
          cycles = [];
          times = [];
          thresholds = ones(1,iNerves)*2000;
          if exist([f,filesep,'cycles.mat'],'file')
            load([fActive,'cycles']);
            dCycles = cycles;
            tCycles = times;
            yCycles = thresholds;
            set(hGraphCycle,'UIContextMenu',hMenuCycle2,'XLim',[1,getCycleMax()]);
            setCameras(hGraphCycle);
          else
            dCycles = [];
            tCycles = [];
            yCycles = ones(1,iNerves)*2000;
            set(hGraphCycle,'UIContextMenu',hMenuCycle);
          end
          dAllNerves = [];
          cbSetTimeScales(hMenuTimeScale,event);
        else
          % TODO (cbSetFolder) warn if neither decimated.mat/spikes.mat/cycles.mat nor nerve*.mat/rawdata*.mat files exist
          warndlg({'This folder has no multichannel data.';'Please choose a new folder.'},'Error reading folder...');
        end
      else
        if exist([f,filesep,'messages.mat'],'file')
          % TODO (cbSetFolder) warn if any data file already exists in the folder, not just messages.mat
          warndlg({'This folder already has multichannel data.';'Please choose a new folder or use the default.'},'Error reading folder...');
        elseif ~strcmp(fActive,[f,filesep])
          if strncmp(getState(),'r',1)
            try
              if ~isempty(dir([fActive,'*']))
                movefile([fActive,'*'],f);
              end
            catch
              lasterror
            end
            try
              rmdir(fActive);
            catch
              lasterror
            end
          else
            set(hToolRecord,'Enable','on');
            set(hToolIdle,'State','on');
            cbSetState(hToolIdle,event);
          end
          fActive = [f,filesep];
          setMessage(['Using folder: ',fActive]);
        end
      end
    end
  end

  % callback for nerve number tools
  function cbSetNerveNumber(source,event)
    if strcmp(get(source,'State'),'on')
      switch source
        case hToolNerve1
          set(hToolNerve2,'State','off');
          iNerves = 1;
          set(hMenuAuxCTH,'Visible','off');
        case hToolNerve2
          set(hToolNerve1,'State','off');
          iNerves = 2;
          if get(hMenuAux,'Value') == 1
            set(hMenuAuxCTH,'Visible','on');
          end
      end
      set(hGraphNerveLine,'Visible','off');
      set(hGraphNerveLine(1:iNerves),'Visible','on');
      set(hGraphCycleLine,'XData',[],'YData',[]);
      hGraphCycleLine = [];
      if strncmp(getState(),'r',1)
        tdtInitializeData(get(hTimer,'UserData'));
      end
    else
      set(source,'State','on');
    end
  end

  % callback for nerve volt scale menu
  function cbSetNerveVoltScales(source,event)
    n = get(source,'UserData');
    set([hGraphNerve,hGraphCycle],'YLim',[-2000,n(get(source,'Value'))]);
    setCameras([hGraphNerve,hGraphCycle]);
    if get(hMenuAux,'Value') == 1
      t = get(hSliderTime,'Value');
      setAux(getTimeMinimum(t),t);
    end
  end

  % callback for print spike frequencies button
  function cbSetOutput(source,event)
    if ~isempty(tSpikes)
      t = get(hSliderTime,'Value');
      tm = getTimeMinimum(t);
      switch get(hMenuOutput,'Value')
        case 1 % spike frequencies
          disp(['Spike Frequencies in Hz from ',getTimeString(tm),' to ',getTimeString(t)]);
          for c = 1:iChannels
            if isempty(sSpikes{c})
              setSpikeSorts(c);
            end
            r = getRange(tSpikes{c},tm,t);
            disp(sprintf('%.3f\t%.3f',sum(sSpikes{c}(r))/(t-tm),length(r)/(t-tm)));
          end
        case 2 % spike times
          disp(sprintf('Spike Times for channel %d in seconds from %s to %s',iChannel,getTimeString(tm),getTimeString(t)));
          if isempty(sSpikes{iChannel})
            setSpikeSorts(iChannel);
          end
          r = getRange(tSpikes{iChannel},tm,t);
          disp(tSpikes{iChannel}(r(sSpikes{iChannel}(r) ~= 0))');
      end
    end
  end

  % callback for changing visible panels
  function cbSetPositions(source,event)
    iMonitors = size(get(0,'MonitorPosition'),1);
    if iMonitors > 1
      set([hToolAnalysis,hToolStatus],'State','on');
    elseif isempty(source) || source == hToolAnalysis
      set(hToolAnalysis,'State','on');
      set(hToolStatus,'State','off');
    elseif source == hToolStatus
      set(hToolAnalysis,'State','off');
      set(hToolStatus,'State','on');
    end
    if strcmp(get(hToolAnalysis,'State'),'on')
      set([hLabelCycle,hTextCycle,hGraphCycle,hLabelThresh,hTextThresh,hGraphThresh,hLabelSort,hGraphSort,hCheckSortEasy,hCheckSortLow, ...
        hCheckSortOutliers,hMenuSortDirection,hGraphAux,hTextAux,hMenuAux,hMenuAuxCTH,hMenuOutput,hButtonOutput],'Visible','on');
      if get(hMenuAux,'Value') > 1 || iNerves < 2
        set(hMenuAuxCTH,'Visible','off');
      end
    else
      set([hLabelCycle,hTextCycle,hGraphCycle,hLabelThresh,hTextThresh,hGraphThresh,hLabelSort,hGraphSort,hCheckSortEasy,hCheckSortLow, ...
        hCheckSortOutliers,hMenuSortDirection,hGraphAux,hTextAux,hMenuAux,hMenuAuxCTH,hMenuOutput,hButtonOutput],'Visible','off');
    end
    if strcmp(get(hToolStatus,'State'),'on')
      set([hBorderRX5,hLabelRX5,hTextRX5Status,hTextRX5Data,hLabelRX5Filter,hMenuRX5FilterHighPass,hLabelRX5FilterDash,hMenuRX5FilterLowPass, ...
        hBorderRP2,hLabelRP2,hTextRP2Status,hTextRP2Data,hBorderUsages,hLabelUsages,hGraphUsages, ...
        hBorderEvents,hLabelEvents,hButtonEvent,hMenuEvents,hTextEvents,hBorderMessages,hLabelMessages,hTextMessages],'Visible','on');
    else
      set([hBorderRX5,hLabelRX5,hTextRX5Status,hTextRX5Data,hLabelRX5Filter,hMenuRX5FilterHighPass,hLabelRX5FilterDash,hMenuRX5FilterLowPass, ...
        hBorderRP2,hLabelRP2,hTextRP2Status,hTextRP2Data,hBorderUsages,hLabelUsages,hGraphUsages, ...
        hBorderEvents,hLabelEvents,hButtonEvent,hMenuEvents,hTextEvents,hBorderMessages,hLabelMessages,hTextMessages],'Visible','off');
    end
    if iMonitors == 1 || isempty(source) || source == hMenuView
      pWindow = get(0,'ScreenSize');
      pSettings  = [0,0,124,180];
      wAControl  = 185;
      wStatus    = 278;
      wStatusRX5 = 154;
      wStatusRP2 = wStatus-wStatusRX5;
      hLabel     =  12;
      if iMonitors == 1
        pWindow = [1,31,pWindow(3),pWindow(4)-83];
        wAnalysis = 400;
        xAControl = pWindow(3)-wAControl-1;
        xAnalysis = xAControl-wAnalysis-1;
        xStatus = pWindow(3)-wStatus-1;
        if strcmp(get(hToolAnalysis,'State'),'on')
          pSettings(1) = pWindow(3)-pSettings(3)-1;
          wSlider = xAnalysis-pSettings(3)-1;
          wGraphs = xAnalysis-3;
        else
          pSettings(1) = xStatus-pSettings(3)-1;
          wSlider = pSettings(1)-1;
          wGraphs = pSettings(1)-3;
        end
      else
        pWindow = [1,0,pWindow(3)*iMonitors,pWindow(4)-52];
        wAnalysis = 620;
        xAnalysis = fix(pWindow(3)*.5);
        xAControl = xAnalysis+wAnalysis+2;
        xStatus = pWindow(3)-wStatus-1;
        pSettings(1) = xAnalysis-pSettings(3)-1;
        wSlider = pSettings(1)-1;
        wGraphs = pSettings(1)-3;
      end
      pSettings(2) = pWindow(4)-pSettings(4)+1;
      hAnGraphs = fix((pSettings(2)-1)/3);
      xStatusRP2 = xStatus+wStatusRX5;
      set(hWindow,               'Position',pWindow);
      % position viewer panel
      set(hBorderSettings,       'Position',pSettings);
      set(hTextClock,            'Position',[pSettings(1)+1  pSettings(2)+153 pSettings(3)-1 18               ]);
      set(hTextState,            'Position',[pSettings(1)+1  pSettings(2)+136 67             18               ]);
      set(hTextStateTime,        'Position',[pSettings(1)+71 pSettings(2)+136 53             18               ]);
      set(hLabelFile,            'Position',[pSettings(1)+22 pSettings(2)+119 58             18               ]);
      set(hTextFile,             'Position',[pSettings(1)+83 pSettings(2)+119 24             18               ]);
      set(hLabelGraphSettings,   'Position',[pSettings(1)+1  pSettings(2)+85  pSettings(3)-1 18               ]);
      set(hLabelTimeScale,       'Position',[pSettings(1)+1  pSettings(2)+63  35             18               ]);
      set(hMenuTimeScale,        'Position',[pSettings(1)+39 pSettings(2)+63  85             22               ]);
      set(hLabelNerveVoltScale,  'Position',[pSettings(1)+1  pSettings(2)+42  35             18               ]);
      set(hMenuNerveVoltScale,   'Position',[pSettings(1)+39 pSettings(2)+42  85             22               ]);
      set(hLabelDataVoltScale,   'Position',[pSettings(1)+1  pSettings(2)+21  35             18               ]);
      set(hMenuDataVoltScale,    'Position',[pSettings(1)+39 pSettings(2)+21  85             22               ]);
      set(hLabelView,            'Position',[pSettings(1)+1  pSettings(2)     35             18               ]);
      set(hMenuView,             'Position',[pSettings(1)+39 pSettings(2)     85             22               ]);
      set(hSliderTime,           'Position',[1               1                wSlider        22               ]);
      set(hMenuTime,             'Position',[wSlider+1       1                pSettings(3)   22               ]);
      set(hLabelNerve,           'Position',[2               pSettings(2)     44             hLabel           ]);
      set(hGraphNerve,           'Position',[1               pSettings(2)     wGraphs        pSettings(4)     ]);
      set(hGraphScope,           'Position',[1               24               wGraphs        pWindow(4)-24    ]);
      % TODO (cbSetPositions) hMenuView - selected channels
      if get(hMenuView,'Value') == 1
        m = 1;
        n = iChannels;
      else
        m = 1+(get(hMenuView,'Value')-2)*16;
        n = m+15;
      end
      set([hLabelData,hGraphData,hGraphDataLine,hGraphSpikes,hGraphSpikesLine],'Visible','off');
      for c = m:n
        y = fix((pSettings(2)-24)*(c-m)/(n-m+1))+24;
        h = fix((pSettings(2)-24)*(c-m+1)/(n-m+1))+24-y;
        set(hLabelData(c),'Position',[2,y,14,min(hLabel,h)]);
        set(hGraphData(c),'Position',[1,y,wGraphs,h]);
        set([hLabelData(c),hGraphData(c),hGraphDataLine(c)],'Visible','on');
        if strcmp(get(hToolStatus,'State'),'on')
          set(hGraphSpikes(c),'Position',[pSettings(1),y,pSettings(3),h]);
          set([hGraphSpikes(c),hGraphSpikesLine(c)],'Visible','on');
        end
      end
      % position analysis panel
      set(hLabelCycle,           'Position',[xAnalysis+2     pSettings(2)     34             hLabel           ]);
      set(hTextCycle,            'Position',[xAnalysis+42    pSettings(2)     wAnalysis-42   hLabel           ]);
      set(hGraphCycle,           'Position',[xAnalysis+1     pSettings(2)     wAnalysis      pSettings(4)     ]);
      set(hLabelThresh,          'Position',[xAnalysis+2     hAnGraphs*2+1    78             hLabel           ]);
      set(hTextThresh,           'Position',[xAnalysis+86    hAnGraphs*2+1    wAnalysis-86   hLabel           ]);
      set(hGraphThresh,          'Position',[xAnalysis+1     hAnGraphs*2+1    wAnalysis      hAnGraphs-1      ]);
      set(hLabelSort,            'Position',[xAnalysis+2     hAnGraphs+1      50             hLabel           ]);
      set(hGraphSort,            'Position',[xAnalysis+1     hAnGraphs+1      wAnalysis      hAnGraphs-1      ]);
      set(hCheckSortEasy,        'Position',[xAControl       hAnGraphs*2-22   wAControl      22               ]);
      set(hCheckSortLow,         'Position',[xAControl       hAnGraphs*2-46   wAControl      22               ]);
      set(hCheckSortOutliers,    'Position',[xAControl       hAnGraphs*2-70   wAControl      22               ]);
      set(hMenuSortDirection,    'Position',[xAControl       hAnGraphs*2-94   wAControl      22               ]);
      set(hGraphAux,             'Position',[xAnalysis+1     1                wAnalysis      hAnGraphs-1      ]);
      set(hTextAux,              'Position',[xAnalysis+2     1                wAnalysis-2    hLabel           ]);
      set(hMenuAux,              'Position',[xAControl       hAnGraphs-22     wAControl      22               ]);
      set(hMenuAuxCTH,           'Position',[xAControl       hAnGraphs-44     wAControl      22               ]);
      set(hMenuOutput,           'Position',[xAControl       25               wAControl      22               ]);
      set(hButtonOutput,         'Position',[xAControl       1                wAControl      22               ]);
      % position status panel
      set(hBorderRX5,            'Position',[xStatus         pSettings(2)     wStatusRX5     pSettings(4)     ]);
      set(hLabelRX5,             'Position',[xStatus+1       pSettings(2)+153 wStatusRX5-1   18               ]);
      set(hTextRX5Status,        'Position',[xStatus+1       pSettings(2)+133 wStatusRX5-1   18               ]);
      set(hTextRX5Data,          'Position',[xStatus+1       pSettings(2)+113 wStatusRX5-1   18               ]);
      set(hLabelRX5Filter,       'Position',[xStatus+1       pSettings(2)+73  wStatusRX5-1   18               ]);
      set(hMenuRX5FilterHighPass,'Position',[xStatus+11      pSettings(2)+49  62             22               ]);
      set(hLabelRX5FilterDash,   'Position',[xStatus+73      pSettings(2)+51  9              18               ]);
      set(hMenuRX5FilterLowPass, 'Position',[xStatus+82      pSettings(2)+49  62             22               ]);
      set(hBorderRP2,            'Position',[xStatusRP2      pSettings(2)     wStatusRP2     pSettings(4)     ]);
      set(hLabelRP2,             'Position',[xStatusRP2+1    pSettings(2)+153 wStatusRP2-1   18               ]);
      set(hTextRP2Status,        'Position',[xStatusRP2+1    pSettings(2)+133 wStatusRP2-1   18               ]);
      set(hTextRP2Data,          'Position',[xStatusRP2+1    pSettings(2)+113 wStatusRP2-1   18               ]);
      set(hBorderUsages,         'Position',[xStatus         pSettings(2)-200 wStatus        199              ]);
      set(hLabelUsages,          'Position',[xStatus+1       pSettings(2)-28  wStatus-1      18               ]);
      set(hGraphUsages,          'Position',[xStatus+25      pSettings(2)-160 wStatus-49     122              ]);
      set(hBorderEvents,         'Position',[xStatus         pSettings(2)-300 wStatus        99               ]);
      set(hLabelEvents,          'Position',[xStatus+1       pSettings(2)-227 wStatus-1      18               ]);
      set(hButtonEvent,          'Position',[xStatus+16      pSettings(2)-249 70             20               ]);
      set(hMenuEvents,           'Position',[xStatus+101     pSettings(2)-249 wStatus-116    20               ]);
      set(hTextEvents,           'Position',[xStatus+16      pSettings(2)-285 wStatus-31     20               ]);
      set(hBorderMessages,       'Position',[xStatus         1                wStatus        pSettings(2)-302 ]);
      set(hLabelMessages,        'Position',[xStatus+1       pSettings(2)-327 wStatus-1      18               ]);
      set(hTextMessages,         'Position',[xStatus+16      17               wStatus-31     pSettings(2)-349 ]);
    end
    setCameras(hWindow);
  end

  % callback for probe number tools
  function cbSetProbeNumber(source,event)
    if strcmp(get(source,'State'),'on')
      switch source
        case hToolProbe1
          set([hToolProbe2,hToolProbe3,hToolProbe4],'State','off');
          set(hMenuView,'Value',1,'String',{'All'});
          iProbes = 1;
        case hToolProbe2
          set([hToolProbe1,hToolProbe3,hToolProbe4],'State','off');
          set(hMenuView,'Value',1,'String',{'All';'1-16';'17-32'});
          iProbes = 2;
        case hToolProbe3
          set([hToolProbe1,hToolProbe2,hToolProbe4],'State','off');
          set(hMenuView,'Value',1,'String',{'All';'1-16';'17-32';'33-48'});
          iProbes = 3;
        case hToolProbe4
          set([hToolProbe1,hToolProbe2,hToolProbe3],'State','off');
          set(hMenuView,'Value',1,'String',{'All';'1-16';'17-32';'33-48';'49-64'});
          iProbes = 4;
      end
      iChannels = iProbes*16;
      if iChannel > iChannels
        setChannel(1);
      end
      cbSetPositions(hMenuView,event);
      if strncmp(getState(),'r',1)
        tdtInitializeData(get(hTimer,'UserData'));
      end
    else
      set(source,'State','on');
    end
  end

  % callback for state tools
  function cbSetState(source,event)
    if strcmp(get(source,'State'),'on')
      state = getState();
      if strncmp(state,'r',1)
        stop(hTimer);
        if source ~= hToolRecord && strcmp(state,'recording')
          setSave(false);
          setMessage(['Timer Period: ',num2str(get(hTimer,'AveragePeriod'))]);
          setMessage(['TDT Buffers: ',num2str(get(hTimer,'UserData'))]);
        end
        if source == hToolIdle
          set([hToolRecord,hToolPlay],'State','off');
          set([hToolFolderOpen,hToolFolderSave,hToolProbe1,hToolProbe2,hToolProbe3,hToolProbe4,hToolNerve1,hToolNerve2],'Enable','on');
          set([hTextState,hTextStateTime],'String','');
          set([hLabelFile,hTextFile],'Visible','off');
          set(hMenuEvents,'Value',1,'String','Other');
          if tdtIsReady()
            tdtSetIdle();
          end
        elseif tdtIsReady()
          set([hToolFolderOpen,hToolFolderSave,hToolProbe1,hToolProbe2,hToolProbe3,hToolProbe4,hToolNerve1,hToolNerve2],'Enable','off');
          set(hTextStateTime,'String','0:00');
          switch source
            case hToolRecord
              set([hToolPlay,hToolIdle],'State','off');
              set(hTextState,'String','Recording:');
              set([hLabelFile,hTextFile],'Visible','on');
              set(hMenuEvents,'String',{'Baseline';'Condition';'Washout';'Other'});
              setFiles(getSaveFileNum('rawdata'),0);
            case hToolPlay
              set([hToolRecord,hToolIdle],'State','off');
              set(hTextState,'String','Playing:');
              set([hLabelFile,hTextFile],'Visible','off');
              set(hMenuEvents,'Value',1,'String','Other');
          end
          dRawData = cell(1,iChannels);
          d = int16(zeros(1,tdtRX5.bufferSize));
          for n = 1:iChannels
            dRawData{n} = d;
          end
          zRawData = 0;
          dDecData = cell(1,iChannels);
          d = int16(zeros(1,2000*tdtRX5.bufferSize/tdtRX5.blockSize));
          for n = 1:iChannels
            dDecData{n} = d;
          end
          iDecData = tdtRX5.blockSize*.5;
          zDecData = 0;
          dSpikes = cell(1,iChannels);
          tSpikes = cell(1,iChannels);
          sSpikes = cell(1,iChannels);
          hSpikes = cell(1,iChannels);
          ySpikes = ones(1,iChannels)*400;
          dNerves = cell(1,iNerves);
          d = int16(zeros(1,tdtRP2.bufferSize));
          for n = 1:iNerves
            dNerves{n} = d;
          end
          zNerves = 0;
          dDecNerves = cell(1,iNerves);
          d = int16(zeros(1,2000*tdtRP2.bufferSize/tdtRP2.blockSize));
          for n = 1:iNerves
            dDecNerves{n} = d;
          end
          iDecNerves = tdtRP2.blockSize*.5;
          dAllNerves = cell(1,iNerves);
          dCycles = cell(1,iNerves);
          tCycles = cell(1,iNerves);
          yCycles = ones(1,iNerves)*2000;
          tdtSetRun();
          setMessage(['State changed to ',get(source,'TooltipString')]);
          tStart = clock;
        else
          set(source,'State','off');
        end
        start(hTimer);
      else
        switch source
          case hToolPlay
            set([hToolRecord,hToolIdle],'State','off');
          case hToolIdle
            set([hToolRecord,hToolPlay],'State','off');
            set(hTextState,'String','Replaying:');
            set([hLabelFile,hTextFile],'Visible','on');
            set(hMenuEvents,'Value',1,'String','Other');
        end
      end
    else
      set(source,'State','on');
    end
  end

  % callback for sorting spikes
  function cbSetSort(source,event)
    t = get(hSliderTime,'Value');
    tm = getTimeMinimum(t);
    setSpikeSorts(iChannel);
    setSpikes(tm,t);
    setAux(tm,t);
  end

  % callback for threshold context menu
  % TODO (cbSetThreshold) improve detection with adaptive thresholding
  % TODO (cbSetThreshold) improve alignment using max of smoothed derivative
  function cbSetThreshold(source,event)
    switch source
      case hMenuCycleSave
        if isempty(dCycles)
          getCycles();
        end
        cycles = dCycles;
        times = tCycles;
        thresholds = yCycles;
        save([fActive,'cycles.mat'],'cycles','times','thresholds');
        set(hGraphCycle,'UIContextMenu',hMenuCycle2);
      case hMenuCycleAuto
        % TODO? (cbSetThreshold) hMenuCycleAuto
        yCycles = ones(1,iNerves)*2000;
        getCycles();
      case hMenuCycleRaise
        for n = 1:iNerves
          yCycles(n) = yCycles(n)+500;
        end
        getCycles();
      case hMenuCycleLower
        for n = 1:iNerves
          yCycles(n) = max(500,yCycles(n)-500);
        end
        getCycles();
      case hMenuCycleEdit
        set(hGraphCycle,'UIContextMenu',hMenuCycle);
      case hMenuThreshSave
        % TODO (cbSetThreshold) save thresholds before trying to pull spikes
        getSpikes();
        spikes = dSpikes;
        times = tSpikes;
        sortcodes = sSpikes;
        thresholds = ySpikes;
        save([fActive,'spikes.mat'],'spikes','times','sortcodes','thresholds');
        set(hGraphThresh,'UIContextMenu',hMenuThresh2);
        % TODO (cbSetThreshold) play a chime to signal the work is finally done!
      case hMenuThreshAutoAll
        for n = find(ySpikes)
          ySpikes(n) = 1.4*getRMS(get(hGraphDataLine(n),'YData'));
        end
      case hMenuThreshAuto
        ySpikes(iChannel) = 1.4*getRMS(get(hGraphThreshLine,'YData'));
      case hMenuThreshRaise
        ySpikes(iChannel) = max(50,ySpikes(iChannel)+20);
      case hMenuThreshLower
        ySpikes(iChannel) = max(50,ySpikes(iChannel)-20);
      case hMenuThreshIgnore
        ySpikes(iChannel) = 0;
      case hMenuThreshEdit
        set(hGraphThresh,'UIContextMenu',hMenuThresh);
    end
    if source == hMenuThreshAutoAll || source == hMenuThreshAuto || source == hMenuThreshRaise || source == hMenuThreshLower || source == hMenuThreshIgnore
      set(hGraphThreshLimitLine,'YData',ones(1,2)*ySpikes(iChannel));
      set(hGraphThreshLimitLineNeg,'YData',-ones(1,2)*ySpikes(iChannel));
    end
  end

  % callback for timeline slider
  function cbSetTime(source,event)
    if ~isstruct(tdtRX5) || ~isstruct(tdtRP2)
      return
    end
    t = get(source,'Value');
    tm = getTimeMinimum(t);
    state = getState();
    d = get(hGraphScope,'UserData');
    dec = false;
    n = t*tdtRX5.samplingRate;
    file = ceil(n/tdtRX5.bufferSize);
    n = n-(file-1)*tdtRX5.bufferSize;
    block = fix(n/tdtRX5.blockSize);
    bmax = tdtRX5.bufferSize/tdtRX5.blockSize;
    if block == 0
      file = file-1;
      block = bmax;
    end
    if d.lastFile > file || (d.lastFile == file && d.lastBlock > block)
      d.lastFile = 1;
      d.lastBlock = 0;
    end
    d.blockUpdate = min(d.blocks,(file-d.lastFile)*bmax+block-d.lastBlock);
    if ~isempty(dDecData) && ~isempty(dDecNerves)
      if (d.dataDecimation == iDecData && d.nerveDecimation == iDecNerves) || (strncmp(state,'r',1) && get(source,'Value') < get(source,'Max'))
        dec = true;
        d.lastFile = file;
        d.lastBlock = block;
        o = fix(tm*2*tdtRX5.samplingRate/tdtRX5.blockSize);
        r = 1:min(2*d.blocks,zDecData-o);
        for c = 1:iChannels
          set(hGraphDataLine(c),'XData',r,'YData',dDecData{c}(r+o));
        end
        o = fix(tm*2*tdtRP2.samplingRate/tdtRP2.blockSize);
        r = 1:min(2*d.blocks,zDecData-o);
        for n = 1:iNerves
          set(hGraphNerveLine(n),'XData',r,'YData',dDecNerves{n}(r+o));
        end
        d.blockIndex = length(r)*.5;
        if d.blockIndex == d.blocks
          d.blockIndex = 0;
        end
      elseif strncmp(state,'r',1)
        size = min([tdtRX5.bufferSize,d.blocks*tdtRX5.blockSize,zDecData*.5*tdtRX5.blockSize]);
        truncate = size < tdtRX5.bufferSize;
        for c = 1:iChannels
          data = circshift(dRawData{c},[0,-zRawData]);
          if truncate
            data = data(1:size);
          end
          data = getPlotDecimation(data,d.dataDecimation);
          set(hGraphDataLine(c),'XData',1:length(data),'YData',data);
        end
        size = min([tdtRP2.bufferSize,d.blocks*tdtRP2.blockSize,zDecData*.5*tdtRP2.blockSize]);
        truncate = size < tdtRP2.bufferSize;
        for c = 1:iNerves
          data = circshift(dNerves{c},[0,-zNerves]);
          if truncate
            data = data(1:size);
          end
          data = getPlotDecimation(data,d.nerveDecimation);
          set(hGraphNerveLine(c),'XData',1:length(data),'YData',data);
        end
        if truncate
          d.blockIndex = 0;
        else
          d.blockIndex = tdtRP2.bufferSize/tdtRP2.blockSize;
        end
      end
      set(hGraphScopeLine,'XData',ones(1,2)*d.blockIndex/d.blocks);
    end
    if strncmp(state,'v',1)
      set(hTextStateTime,'String',getTimeString(t));
      if dec
        setFiles(d.lastFile,0);
      else
        n = 0;
        files = [];
        blocks = cell(0);
        while n < d.blocks && (d.lastFile < file || (d.lastFile == file && d.lastBlock < block))
          b = [];
          while n < d.blocks && block > 0 && (d.lastFile < file || (d.lastFile == file && d.lastBlock < block))
            n = n+1;
            b = [block,b];
            block = block-1;
          end
          files = [file,files];
          blocks = {b,blocks{:}};
          if block == 0
            file = file-1;
            block = bmax;
          end
        end
        if n == d.blocks
          d.blockIndex = 0;
        end
        for n = 1:length(files)
          b = length(blocks{n});
          for c = 1:iChannels
            setData(hGraphDataLine(c),getPlotDecimation(getRawData(c,files(n),blocks{n}),d.dataDecimation),d.blockIndex,d.blocks,b);
          end
          for c = 1:iNerves
            setData(hGraphNerveLine(c),getPlotDecimation(getNerveData(c,files(n),blocks{n}),d.nerveDecimation),d.blockIndex,d.blocks,b);
          end
          d.blockIndex = d.blockIndex+b;
          if d.blockIndex >= d.blocks
            d.blockIndex = d.blockIndex-d.blocks;
          end
          set(hGraphScopeLine,'XData',ones(1,2)*d.blockIndex/d.blocks);
          setDraw();
          if n == length(files)
            d.lastFile = files(n);
            d.lastBlock = blocks{n}(b);
          end
        end
      end
    end
    set(hGraphScope,'UserData',d);
    set(hGraphThreshLine,'XData',get(hGraphDataLine(iChannel),'XData'),'YData',get(hGraphDataLine(iChannel),'YData'));
    setCycles(tm,t);
    setSpikes(tm,t);
    setAux(tm,t);
  end

  % callback for time scale menu
  function cbSetTimeScales(source,event)
    set([hGraphNerveLine,hGraphDataLine,hGraphThreshLine],'XData',[],'YData',[]);
    set(hGraphScopeLine,'XData',[0,0]);
    if ~isstruct(tdtRX5) || ~isstruct(tdtRP2)
      return
    end
    t = get(source,'UserData');
    t = t(get(source,'Value'));
    b = tdtRX5.blockSize/tdtRX5.samplingRate;
    bmax = tdtRX5.bufferSize/tdtRX5.blockSize;
    if t == 0
      set(hSliderTime,'Visible','off','Value',get(hSliderTime,'Max'));
      if isempty(dDecData)
        iDecData = tdtRX5.blockSize*.5;
        iDecNerves = tdtRP2.blockSize*.5;
        m = min(length(dir([fActive,'rawdata*.mat'])),length(dir([fActive,'nerve1*.mat'])));
        set([hGraphData,hGraphNerve],'XLim',[1,m*tdtRP2.bufferSize/iDecNerves]);
        setCameras([hGraphData,hGraphNerve]);
        for f = 1:m
          setFiles(f,1);
          for n = 1:iChannels
            setData(hGraphDataLine(n),getPlotDecimation(dRawData{n},iDecData),f-1,m,1);
          end
          for n = 1:iNerves
            setData(hGraphNerveLine(n),getPlotDecimation(dNerves{n},iDecNerves),f-1,m,1);
          end
          setDraw();
        end
        dDecData = cell(1,iChannels);
        for n = 1:iChannels
          dDecData{n} = get(hGraphDataLine(n),'YData');
        end
        dDecNerves = cell(1,iNerves);
        for n = 1:iNerves
          dDecNerves{n} = get(hGraphNerveLine(n),'YData');
        end
        decimated = struct;
        decimated.dataDecimation = iDecData;
        decimated.nerveDecimation = iDecNerves;
        decimated.data = dDecData;
        decimated.nerves = dDecNerves;
        save([fActive,'decimated.mat'],'decimated');
        zDecData = length(dDecData{1});
        m = zDecData*.5*b;
        set(hSliderTime,'Value',m,'Max',m);
      end
      set([hGraphData,hGraphNerve,hGraphThresh],'XLim',[1,zDecData]);
      s = zDecData*.5;
    else
      s = ceil(t/b);
    end
    d = struct;
    d.lastFile = 1;
    d.lastBlock = 0;
    d.blocks = s;
    d.blockIndex = 0;
    d.blockUpdate = 0;
    d.dataDecimation = 1;
    d.nerveDecimation = 1;
    if t == 0
      f = ceil(s/bmax);
      d.lastFile = f;
      d.lastBlock = s-bmax*(f-1);
      d.dataDecimation = iDecData;
      d.nerveDecimation = iDecNerves;
    else
      state = getState();
      nData = s*tdtRX5.blockSize;
      if strncmp(state,'v',1) && ~isempty(dDecData) && isempty(dir([fActive,'rawdata*.mat']))
        d.dataDecimation = iDecData;
      else
        d.dataDecimation = getMaxFactor(tdtRX5.blockSize,nData/1152);
        if mod(tdtRX5.blockSize/d.dataDecimation,2) == 1
          d.dataDecimation = getMaxFactor(tdtRX5.blockSize*.5,nData/1152);
        end
      end
      set([hGraphData,hGraphThresh],'XLim',[1,nData/d.dataDecimation]);
      nNerve = s*tdtRP2.blockSize;
      if strncmp(state,'v',1) && ~isempty(dDecNerves) && isempty(dir([fActive,'nerve*.mat']))
        d.nerveDecimation = iDecNerves;
      else
        d.nerveDecimation = getMaxFactor(tdtRP2.blockSize,nNerve/1152);
        if mod(tdtRP2.blockSize/d.nerveDecimation,2) == 1
          d.nerveDecimation = getMaxFactor(tdtRP2.blockSize*.5,nNerve/1152);
        end
      end
      set(hGraphNerve,'XLim',[1,nNerve/d.nerveDecimation]);
      t = s*b;
      m = get(hSliderTime,'Max');
      if t >= m
        set(hSliderTime,'Visible','off','Value',m);
      else
        set(hSliderTime,'Visible','on','Value',max(get(hSliderTime,'Value'),t),'Min',t,'SliderStep',[min(.1,b/(m-t)),t/(m-t)]);
      end
    end
    set(hGraphScope,'UserData',d);
    set([hGraphThreshLimitLine,hGraphThreshLimitLineNeg],'XData',get(hGraphThresh,'XLim'));
    setCameras([hGraphData,hGraphNerve,hGraphThresh]);
    cbSetTime(hSliderTime,event);
  end

  % callback for main update timer
  function cbTimer(source,event)
    persistent tUpdate;
    set(hTextClock,'String',datestr(clock));
    state = getState();
    if strncmp(state,'r',1)
      tUpdate = tUpdate+get(source,'InstantPeriod');
      if isempty(tUpdate)
        tUpdate = 0;
      end
      if ~strcmp(state,'recording,offline')
        if ~strcmp(state,'recording,idle')
          set(hTextStateTime,'String',getTimeString(etime(clock,tStart)));
          rRawData = 2*tdtGetSync(tdtRX5);
          rNerves  = 2*tdtGetSync(tdtRP2);
          if (rRawData >= zRawData+tdtRX5.blockSize || rRawData < zRawData) && (rNerves >= zNerves+tdtRP2.blockSize || rNerves < zNerves)
            rx5block = zRawData+(1:tdtRX5.blockSize);
            rp2block = zNerves+(1:tdtRP2.blockSize);
            for n = 1:iChannels
              dRawData{n}(rx5block) = tdtGetData(tdtRX5,[tdtRX5.data,'~',int2str(n)],zRawData);
              dDecData{n}(zDecData+(1:2)) = getPlotDecimation(dRawData{n}(rx5block),iDecData);
            end
            for n = 1:iNerves
              dNerves{n}(rp2block) = tdtGetData(tdtRP2,[tdtRP2.data(1:4),int2str(n)],zNerves);
              dDecNerves{n}(zDecData+(1:2)) = getPlotDecimation(dNerves{n}(rp2block),iDecNerves);
              % TODO (cbTimer) find cycles online
%              dAllNerves{n}(end+(1:tdtRP2.blockSize)) = dNerves{n}(rp2block);
            end
            zDecData = zDecData+2;
            if zDecData == length(dDecData{1})
              d = int16(zeros(1,2000*tdtRX5.bufferSize/tdtRX5.blockSize));
              for n = 1:iChannels
                dDecData{n} = [dDecData{n},d];
              end
              d = int16(zeros(1,2000*tdtRP2.bufferSize/tdtRP2.blockSize));
              for n = 1:iNerves
                dDecNerves{n} = [dDecNerves{n},d];
              end
            end
            d = get(hGraphScope,'UserData');
            b = tdtRX5.blockSize/tdtRX5.samplingRate;
            t = d.blocks*b;
            m = zDecData*.5*b;
            % TODO (cbTimer) handle showing data when time scale == all data
            if get(hSliderTime,'Value') == get(hSliderTime,'Max')
              for n = 1:iChannels
                setData(hGraphDataLine(n),getPlotDecimation(dRawData{n}(rx5block),d.dataDecimation),d.blockIndex,d.blocks,1);
              end
              set(hGraphThreshLine,'XData',get(hGraphDataLine(iChannel),'XData'),'YData',get(hGraphDataLine(iChannel),'YData'));
              for n = 1:iNerves
                setData(hGraphNerveLine(n),getPlotDecimation(dNerves{n}(rp2block),d.nerveDecimation),d.blockIndex,d.blocks,1);
              end
              d.blockIndex = d.blockIndex+1;
              if d.blockIndex >= d.blocks
                d.blockIndex = d.blockIndex-d.blocks;
              end
              d.blockUpdate = 1;
              set(hGraphScopeLine,'XData',ones(1,2)*d.blockIndex/d.blocks);
              set(hGraphScope,'UserData',d);
              set(hSliderTime,'Value',m,'Max',m);
            else
              set(hSliderTime,'Max',m);
            end
            if t >= m
              set(hSliderTime,'Visible','off');
            else
              set(hSliderTime,'Visible','on','Value',max(get(hSliderTime,'Value'),t),'Min',t,'SliderStep',[min(.1,b/(m-t)),t/(m-t)]);
            end
            zRawData = zRawData+tdtRX5.blockSize;
            zNerves = zNerves+tdtRP2.blockSize;
            if zRawData >= length(dRawData{1}) && zNerves >= length(dNerves{1})
              zRawData = zRawData-length(dRawData{1});
              zNerves = zNerves-length(dNerves{1});
              if strcmp(state,'recording')
                setSave(true);
              else
                setMessage('Filled buffers');
              end
            end
          end
          % TODO (cbTimer) get spikes and cycles from recent data
          % TODO (cbTimer) show spikes on spike and sort graphs
          % TODO (cbTimer) update auxiliary graph
        end
        if tUpdate >= 1
          tUpdate = 0;
          set(hTextRX5Status,'String',['RX5 is ',tdtGetStatus(tdtRX5.device),'.']);
          set(hTextRP2Status,'String',['RP2 is ',tdtGetStatus(tdtRP2.device),'.']);
          u = tdtGetUsage();
          set(hGraphUsageRX5M,'YData',[0,0,u(1),u(1)],'FaceVertexCData',[.4,1,.4;.4,1,.4;.4+.006*u(1),1-.006*u(1),.4;.4+.006*u(1),1-.006*u(1),.4]);
          set(hGraphUsageRX5A,'YData',[0,0,u(2),u(2)],'FaceVertexCData',[.4,1,.4;.4,1,.4;.4+.006*u(2),1-.006*u(2),.4;.4+.006*u(2),1-.006*u(2),.4]);
          set(hGraphUsageRP2, 'YData',[0,0,u(3),u(3)],'FaceVertexCData',[.4,1,.4;.4,1,.4;.4+.006*u(3),1-.006*u(3),.4;.4+.006*u(3),1-.006*u(3),.4]);
        end
      else
        if strcmp(get(hToolIdle,'State'),'off')
          set(hToolIdle,'State','on');
          cbSetState(hToolIdle,event);
        end
        set(hTextRX5Status,'String','RX5 is disconnected.');
        set(hTextRP2Status,'String','RP2 is disconnected.');
        set([hGraphUsageRX5M,hGraphUsageRX5A,hGraphUsageRP2],'YData',[0,0,0,0]);
        if tUpdate >= 5
          tUpdate = 0;
          setTDTDefaults(0);
        end
      end
    elseif strcmp(state,'viewing')
      set(hSliderTime,'Value',min(get(hSliderTime,'Max'),get(hSliderTime,'Value')+get(source,'InstantPeriod')*2^(get(hMenuTime,'Value')-4)));
      cbSetTime(hSliderTime,event);
    end
  end

  % returns a color from list
  function c = getColor(n)
    colors = {'white','red','green','blue','magenta','cyan','yellow'};
    c = colors{1+mod(n-1,length(colors))};
  end

  % returns max length of cycle
  function x = getCycleMax()
    x = 2;
    for n = 1:iNerves
      for c = 1:length(dCycles{n})
        x = max(x,length(dCycles{n}{c}));
      end
    end
  end

  % updates list of cycles
  % TODO (getCycles) use area threshold instead of 0.5*gap
  % TODO (getCycles) let user change gap from default of 5 sec (1905 points)
  function getCycles()
    if ~isstruct(tdtRP2)
      return
    end
    if isempty(dAllNerves)
      dAllNerves = cell(1,iNerves);
      s = str2double(get(hTextFile,'String'));
      for n = 1:iNerves
        nerve = [];
        files = dir(sprintf('%snerve%d*.mat',fActive,n));
        for f = 1:length(files)
          setFiles(f,0);
          setDraw();
          load([fActive,files(f).name]);
          dAllNerves{n}(end+(1:length(nerve))) = nerve;
        end
        for w = find(dAllNerves{n} < -2000)
          dAllNerves{n}(w) = dAllNerves{n}(w)+65535;
        end
      end
      setFiles(s,0);
    end
    dCycles = cell(1,iNerves);
    tCycles = cell(1,iNerves);
    gap = 1905;
    halfgap = fix(.5*gap);
    for n = 1:iNerves
      x = 1;
      y = find(dAllNerves{n}(1:end-gap) > yCycles(n));
      while x+halfgap < length(y)
        if y(x+halfgap)-y(x) < gap
          t = y(x);
          x = x+halfgap;
          while x+1 < length(y) && y(x+1)-y(x) < gap
            x = x+1;
          end
          dCycles{n}{end+1} = dAllNerves{n}(max(1,t-gap):min(length(dAllNerves{n}),y(x)+gap));
          tCycles{n}(end+1) = t/tdtRP2.samplingRate;
        end
        x = x+1;
      end
    end
    set(hGraphCycle,'XLim',[1,getCycleMax()]);
    setCameras(hGraphCycle);
    t = get(hSliderTime,'Value');
    tm = getTimeMinimum(t);
    setCycles(tm,t);
    setAux(tm,t);
  end

  % returns a list of factors for a given number
  function factors = getFactors(n)
    factors = unique(factor(n));
    if length(factors) > 1
      for f = factors
        factors = unique([factors,n/f,getFactors(n/f)]);
      end
    end
  end

  % returns icon
  function cdata = getIcon(s)
    cdata = [];
    if exist(fullfile('Icons',[s '.png']),'file') == 2
      [cdata,map,alpha] = imread(fullfile('Icons',[s '.png']));
    end
    if isempty(cdata)
      return;
    end
    if isempty(map)
      cdata = double(cdata);
      cdata = cdata/255;
    else
      cdata = ind2rgb(cdata,map);
    end
    if isempty(alpha)
      alpha = ~(cdata(:,:,1) == cdata(1,1,1) & cdata(:,:,2) == cdata(1,1,2) & cdata(:,:,3) == cdata(1,1,3));
    end
    r = cdata(:,:,1);
    r(alpha == 0) = NaN;
    g = cdata(:,:,2);
    g(alpha == 0) = NaN;
    b = cdata(:,:,3);
    b(alpha == 0) = NaN;
    cdata = cat(3,r,g,b);
  end

  % returns greatest factor still less than a given maximum
  function m = getMaxFactor(n,x)
    f = getFactors(n);
    m = max([1,f(f <= x)]);
  end

  % returns nerve data from file
  function data = getNerveData(channel,file,blocks)
    data = [];
    if isstruct(tdtRP2)
      setFiles(file,1);
      try
        data = dNerves{channel}(((blocks(1)-1)*tdtRP2.blockSize)+(1:length(blocks)*tdtRP2.blockSize));
      catch
        lasterror
        data = dNerves{channel}((blocks(1)-1)*tdtRP2.blockSize+1:length(dNerves{channel}));
        data = [data,zeros(1,length(blocks)*tdtRP2.blockSize-length(data))];
      end
    end
  end

  % returns plot decimation of data
  function p = getPlotDecimation(data,decimation)
    if decimation > 1
      lengthp = length(data)/decimation;
      p = zeros(1,lengthp);
      j = 1:2:lengthp;
      k = 2:2:lengthp;
      l = 1:2*decimation:length(data);
      m = l+(2*decimation-1);
      for n = 1:length(k)
        p(j(n)) = min(data(l(n):m(n)));
        p(k(n)) = max(data(l(n):m(n)));
      end
    else
      p = data;
    end
  end

  % returns index of data near given value, assumes ordered list
  function p = getPoint(data,n)
    l = 1;
    h = length(data);
    p = fix((l+h)*.5);
    while l <= h
      if data(p) < n
        l = p+1;
      elseif n < data(p)
        h = p-1;
      else
        break
      end
      p = fix((l+h)*.5);
    end
    if p == 0
      p = 1;
    end
  end

  % returns indices of data within range, assumes ordered list
  function r = getRange(data,low,high)
    r = [];
    if ~isempty(data)
      l = getPoint(data,low);
      if data(l) < low
        l = l+1;
      end
      h = getPoint(data,high);
      if data(h) > high
        h = h-1;
      end
      if l <= h
        r = l:h;
      end
    end
  end

  % returns raw multichannel data from file
  function data = getRawData(channel,file,blocks)
    data = [];
    if isstruct(tdtRX5)
      setFiles(file,1);
      try
        data = dRawData{channel}(((blocks(1)-1)*tdtRX5.blockSize)+(1:length(blocks)*tdtRX5.blockSize));
      catch
        lasterror
        data = dRawData{channel}((blocks(1)-1)*tdtRX5.blockSize+1:length(dRawData{channel}));
        data = [data,zeros(1,length(blocks)*tdtRX5.blockSize-length(data))];
      end
    end
  end

  % returns root-mean-square (RMS) value of data
  function rms = getRMS(data)
    rms = length(data);
    if rms ~= 0
      rms = sqrt(sum(data.*data)/rms);
    end
  end

  % returns next number of file
  function n = getSaveFileNum(s)
    n = length(dir([fActive,s,'*.mat']));
    if n == 0
      n = 1;
    end
    while exist(sprintf('%s%s-%03d.mat',fActive,s,n),'file')
      n = n + 1;
    end
  end

  % updates list of spikes
  function getSpikes()
    if ~isstruct(tdtRX5)
      return
    end
    dSpikes = cell(1,iChannels);
    tSpikes = cell(1,iChannels);
    sSpikes = cell(1,iChannels);
    hSpikes = cell(1,iChannels);
    rawdata = [];
    files = dir([fActive,'rawdata*.mat']);
    s = str2double(get(hTextFile,'String'));
    o = 0;
    t = 1/tdtRX5.samplingRate;
    for f = 1:length(files)
      setFiles(f,0);
      setDraw();
      load([fActive,files(f).name]);
      for c = 1:iChannels
        times = getSpikeTimes(rawdata{c},ySpikes(c),62);
        if ~isempty(times)
          dSpikes{c}(:,end+(1:length(times))) = getSpikeSnips(rawdata{c},times,62);
          tSpikes{c}(end+(1:length(times))) = (times+o)*t;
        end
      end
      o = o+tdtRX5.bufferSize;
    end
    % TODO (getSpikes) fill sSpikes with sort codes, set saved status
    setFiles(s,0);
    t = get(hSliderTime,'Value');
    tm = getTimeMinimum(t);
    setSpikes(tm,t);
    setAux(tm,t);
  end

  % picks out spike shape from spike snip: time and height
  function [t,h] = getSpikeShape(spike)
    % find threshold index
    n = fix(.5*length(spike));
    % invert spike if negative threshold
    if spike(n) < 0
      s = -spike;
    else
      s = spike;
    end
    t = [];
    e = [];
    % find previous extreme
    t(1) = getPreviousStart(n);
    n = getStillNegative(t(1));
    if t(1) > 0
      [e(1),d] = min(s(t(1):n));
      t(1) = t(1)+d-1;
    end
    % find threshold extreme
    t(2) = n+1;
    n = getStillPositive(t(2));
    [e(2),d] = max(s(t(2):n));
    t(2) = t(2)+d-1;
    % find next extreme
    if n < length(s)
      t(3) = n+1;
      n = getStillNegative(t(3));
      [e(3),d] = min(s(t(3):n));
      t(3) = t(3)+d-1;
    end
    % find next extreme
    if n < length(s)
      t(4) = n+1;
      n = getStillPositive(t(4));
      [e(4),d] = max(s(t(4):n));
      t(4) = t(4)+d-1;
    end
    % find max height between extremes
    d = diff(e);
    [h,n] = max(abs(d));
    h = d(n);
    % find time between those extremes
    t = t(n+1)-t(n);

    function x = getPreviousStart(x)
      while x > 0 && s(x) > 0
        x = x-1;
      end
      while x > 1 && s(x-1) < 0
        x = x-1;
      end
    end

    function x = getStillNegative(x)
      while x < length(s) && s(x+1) < 0
        x = x+1;
      end
    end

    function x = getStillPositive(x)
      while x < length(s) && s(x+1) > 0
        x = x+1;
      end
    end

  end

  % returns snips out of data
  function snips = getSpikeSnips(data,times,width)
    snips = [];
    if ~isempty(times)
      halfwidth = fix(width*.5);
      for n = 1:length(times)
        snips(:,end+1) = (data(times(n)-halfwidth+(1:width)))';
      end
    end
  end

  % returns times when threshold is exceeded
  function times = getSpikeTimes(data,thresh,width)
    times = [];
    if thresh ~= 0
      halfwidth = fix(width*.5);
      n = halfwidth;
      while n < length(data)-halfwidth
        if thresh < abs(data(n))
          times(end+1) = n;
          n = n+halfwidth;
        else
          n = n+1;
        end
      end
    end
  end

  % returns state (recording:preview,idle,offline;viewing:idle)
  function state = getState()
    try
      if strcmp(get(hToolRecord,'State'),'on')
        state = 'recording';
      elseif strcmp(get(hToolRecord,'Enable'),'on')
        if strcmp(get(hToolPlay,'State'),'on')
          state = 'recording,preview';
        elseif tdtIsReady()
          state = 'recording,idle';
        else
          state = 'recording,offline';
        end
      elseif strcmp(get(hToolPlay,'State'),'on')
        state = 'viewing';
      else
        state = 'viewing,idle';
      end
    catch
      lasterror
      state = [];
    end
  end

  % returns time minimum of visible data
  function tm = getTimeMinimum(t)
    tm = get(hMenuTimeScale,'UserData');
    tm = tm(get(hMenuTimeScale,'Value'));
    if tm ~= 0
      tm = max(0,t-tm);
    end
  end

  % returns time string hh:mm:ss
  function s = getTimeString(t)
    h = fix(t/3600);
    tm = t-h*3600;
    m = fix(tm/60);
    if h > 0
      s = sprintf('%d:%02d:%02d',h,m,fix(tm-m*60));
    else
      s = sprintf('%d:%02d',m,fix(tm-m*60));
    end
  end

  % updates auxiliary graph
  function setAux(tm,t)
    set(hTextAux,'String','');
    delete(findobj(hGraphAux,'Type','line'));
    delete(findobj(hGraphAux,'Type','patch'));
    set(hWindow,'CurrentAxes',hGraphAux);
    switch get(hMenuAux,'Value')
      case 1 % Cycle Triggered Histogram
        % TODO (setAux) cth controls - bin number, count scale, set time rather than cycle length, mouseover information
        % TODO (setAux) cth analysis - define pre,E,I,post boundaries, classify, display eta^2 statistic of respiratory relation
        % TODO (setAux) check for sSpikes sort codes
        if ~isempty(tCycles) && ~isempty(tSpikes) && isstruct(tdtRP2)
          bins = 40;
          times = [];
          cycleLength = 2;
          n = 1;
          if iNerves > 1
            n = get(hMenuAuxCTH,'Value');
          end
          cycleIndices = getRange(tCycles{n},tm,t);
          if ~isempty(cycleIndices)
            times = cell(1,length(cycleIndices));
            cycleAverage = zeros(1,getCycleMax(),'int32');
            for i = cycleIndices
              cycle = dCycles{n}(i);
              for x = 1:length(cycle{1})
                cycleAverage(x) = cycleAverage(x)+int32(cycle{1}(x));
              end
              cycleTime = tCycles{n}(i);
              times{i-cycleIndices(1)+1} = tSpikes{iChannel}(getRange(tSpikes{iChannel},cycleTime-5,cycleTime-5+(length(cycle{1})/tdtRP2.samplingRate)))-cycleTime;
            end
            cycleAverage = cycleAverage/length(cycleIndices);
            cycleAverage = cycleAverage(1:find(cycleAverage,1,'last'));
            line(1:length(cycleAverage),cycleAverage,'Color',getColor(n));
            cycleLength = max(cycleLength,length(cycleAverage));
          end
          ylim = get(hGraphCycle,'YLim');
          set(hGraphAux,'XLim',[1,cycleLength],'YLim',[0,ylim(2)]);
          cycleLength = cycleLength/tdtRP2.samplingRate;
          s = cycleLength*.5/(bins-1);
          x = linspace(-5+s,cycleLength-5-s,bins);
          if ~isempty(times)
            h = zeros(length(times),bins);
            for i = 1:length(times)
              h(i,:) = hist(times{i},x);
            end
            [p,t,s] = anova1(h,[],'off');
            eta = t{2,3}*t{2,5}/(t{2,3}*t{2,5}+t{3,3}); % dfb*F/(dfb*F+dfw), equation 2 in Orem and Dick, 1983
            cth = s.means;
            if max(cth) > 0
              x = (x+5)*tdtRP2.samplingRate;
              s = .9*x(1);
              cth = cth*(ylim(2)/max(cth));
              for n = find(cth)
                patch([x(n)-s,x(n)+s,x(n)+s,x(n)-s],[0,0,cth(n),cth(n)],getColor(iChannel));
              end
              set(hTextAux,'String',sprintf('eta^2 = %.2f',eta));
            end
          end
        end
      case 2 % Spike Triggered Histogram
        % TODO (setAux) spike triggered histogram
      case 3 % Interspike Intervals
        % TODO (setAux) isi controls - time scale, normal/log?, mouseover information
        % TODO (setAux) check for sSpikes sort codes
        if ~isempty(tSpikes)
          isi = diff(tSpikes{iChannel}(getRange(tSpikes{iChannel},tm,t)));
          if ~isempty(isi)
            line(1:length(isi),isi,'Color',getColor(iChannel));
            set(hGraphAux,'XLim',[1,length(isi)],'YLim',[0,max(isi)]);
          end
        end
      case 4 % Interspike Intervals Histogram
        % TODO (setAux) isih controls - bin number, time and count scales, normal/log?, mouseover information
        % TODO (setAux) isih analysis - burster/tonic, burst rate, spike rate in bursts, spikes per burst
        % TODO (setAux) check for sSpikes sort codes
        if ~isempty(tSpikes)
          bins = 25;
          s = 3/(bins-1);
          [h,x] = hist(log10(diff(tSpikes{iChannel}(getRange(tSpikes{iChannel},tm,t)))),linspace(-3,3,bins));
          if max(h) > 0
            for n = find(h)
              h(n) = log10(h(n))+1;
              patch([x(n)-s,x(n)+s,x(n)+s,x(n)-s],[0,0,h(n),h(n)],getColor(iChannel));
            end
            set(hGraphAux,'XLim',[-3,3],'YLim',[0,max(h)]);
          end
        end
      case 5 % Firing Rate Histogram
        % TODO (setAux) firing rate histogram
    end
    setCameras(hGraphAux);
  end

  % updates camera positions for graphs
  function setCameras(h)
    if strcmp(get(hWindow,'Visible'),'on')
      set(findobj(h,'Type','axes','-and','Visible','on'),'CameraPositionMode','auto',  'CameraTargetMode','auto',  'CameraUpVectorMode','auto');
      setDraw();
      set(findobj(h,'Type','axes','-and','Visible','on'),'CameraPositionMode','manual','CameraTargetMode','manual','CameraUpVectorMode','manual');
    end
  end

  % sets channel for analysis panel
  % TODO (setChannel) only update spike sort, not all of the single spike graphs
  function setChannel(n)
    if isempty(iChannel) || iChannel ~= n
      set(hLabelData(iChannel),'FontWeight','default');
      set(hLabelData(n),       'FontWeight','bold');
      set([hGraphData(iChannel),hGraphSpikes(iChannel)],'Color','default');
      set([hGraphData(n),       hGraphSpikes(n)],       'Color',[.2,.2,.2]);
      set([hGraphThreshLine,hGraphSortLine],'XData',[],'YData',[],'Color',getColor(n));
      set(hGraphThreshLine,'XData',get(hGraphDataLine(n),'XData'),'YData',get(hGraphDataLine(n),'YData'));
      set([hGraphThreshLimitLine,hGraphThreshLimitLineNeg],'XData',get(hGraphThresh,'XLim'));
      set(hGraphThreshLimitLine,   'YData',ones(1,2)*ySpikes(n));
      set(hGraphThreshLimitLineNeg,'YData',-ones(1,2)*ySpikes(n));
      iChannel = n;
      v = get(hMenuView,'Value');
      if v > 1
        m = 1+(v-2)*16;
        if n < m
          set(hMenuView,'Value',v-1);
          cbSetPositions(hMenuView,[]);
        elseif n > m+15
          set(hMenuView,'Value',v+1);
          cbSetPositions(hMenuView,[]);
        end
      end
      t = get(hSliderTime,'Value');
      tm = getTimeMinimum(t);
      setSpikes(tm,t);
      setAux(tm,t);
    end
  end

  % sets cycles of nerve activity
  function setCycles(tm,t)
    set(hTextCycle,'String','');
    set(hGraphCycleLine,'XData',[],'YData',[]);
    hGraphCycleLine = [];
    if ~isempty(dCycles)
      set(hWindow,'CurrentAxes',hGraphCycle);
      s = 'Average Interval';
      for n = 1:iNerves
        r = getRange(tCycles{n},tm,t);
        for x = r
          hGraphCycleLine(end+1) = line('Color',getColor(n),'XData',1:length(dCycles{n}{x}),'YData',dCycles{n}{x});
        end
        if length(r) > 1
          abi = (tCycles{n}(length(r))-tCycles{n}(1))/(length(r)-1);
          if abi > 99
            s = sprintf('%s, Nerve %d: %.1f min (%.0f)',s,n,abi/60,length(r));
          else
            s = sprintf('%s, Nerve %d: %.0f sec (%.0f)',s,n,abi,length(r));
          end
          set(hTextCycle,'String',s);
        end
      end
    end
  end

  % adds data to a line
  function setData(line,data,offset,total,num)
    d = get(line,'YData');
    b = length(data)/num;
    d(offset*b+(1:length(data))) = data;
    if offset+num > total
      x = total*b;
      d = [d(x+1:end),d(end-x:x)];
    end
    set(line,'XData',1:length(d),'YData',d);
  end

  % attempts to draw, catches any errors
  function setDraw()
    try
      drawnow();
    catch
      lasterror
    end
  end

  % sets buffers to data from files
  function setFiles(file,open)
    set(hTextFile,'String',int2str(file));
    if open && file ~= get(hTextFile,'UserData')
      set(hTextFile,'UserData',file);
      dRawData = cell(1,iChannels);
      rawdata = [];
      r = dir([fActive,'rawdata*.mat']);
      if length(r) >= file
        load([fActive,r(file).name]);
        dRawData = rawdata;
      end
      dNerves = cell(1,iNerves);
      nerve = [];
      n1 = dir([fActive,'nerve1*.mat']);
      if length(n1) >= file
        load([fActive,n1(file).name]);
        dNerves{1} = nerve;
      end
      n2 = dir([fActive,'nerve2*.mat']);
      if length(n2) >= file
        load([fActive,n2(file).name]);
        dNerves{2} = nerve;
      end
    end
  end

  % sets GUI defaults
  function setGUIDefaults()
    set(0, ...
      'DefaultFigureDockControls',           'off', ...
      'DefaultFigureMenuBar',                'none', ...
      'DefaultFigureNumberTitle',            'off', ...
      'DefaultFigureRenderer',               'OpenGL', ...
      'DefaultFigureRendererMode',           'manual', ...
      'DefaultFigureResize',                 'off', ...
      'DefaultFigureVisible',                'off', ...
      'DefaultUIControlBackgroundColor',     'black', ...
      'DefaultUIControlForegroundColor',     'white', ...
      'DefaultUIControlFontSize',            8, ...
      'DefaultUIControlHorizontalAlignment', 'left', ...
      'DefaultUIControlValue',               1, ...
      'DefaultAxesALimMode',                 'manual', ...
      'DefaultAxesCLimMode',                 'manual', ...
      'DefaultAxesColor',                    'black', ...
      'DefaultAxesTickDirMode',              'manual', ...
      'DefaultAxesUnits',                    'pixels', ...
      'DefaultAxesXColor',                   [.4,.4,.4], ...
      'DefaultAxesXLimMode',                 'manual', ...
      'DefaultAxesXTick',                    [], ...
      'DefaultAxesXTickLabelMode',           'manual', ...
      'DefaultAxesXTickMode',                'manual', ...
      'DefaultAxesYColor',                   [.4,.4,.4], ...
      'DefaultAxesYLimMode',                 'manual', ...
      'DefaultAxesYTick',                    [], ...
      'DefaultAxesYTickLabelMode',           'manual', ...
      'DefaultAxesYTickMode',                'manual', ...
      'DefaultAxesZColor',                   [.4,.4,.4], ...
      'DefaultAxesZLimMode',                 'manual', ...
      'DefaultAxesZTick',                    [], ...
      'DefaultAxesZTickLabelMode',           'manual', ...
      'DefaultAxesZTickMode',                'manual', ...
      'DefaultLineColor',                    'white', ...
      'DefaultLineEraseMode',                'normal', ...
      'DefaultLineHitTest',                  'off', ...
      'DefaultLineXData',                    [], ...
      'DefaultLineYData',                    []);
  end

  % sets messages text in status panel
  function setMessage(s)
    state = getState();
    if strncmp(state,'r',1)
      m = get(hTextMessages,'UserData');
    else
      m = get(hTextMessages,'String');
    end
    if isempty(m) || strcmp(s,'Program started')
      m = [get(hTextClock,'String'),'  ',s];
    elseif ~ishandle(hWindow) || strncmp(s,'Program',7) || strncmp(s,'State',5) || strncmp(s,'Using',5) || (~strcmp(state,'recording') && ~strcmp(state,'recording,preview'))
      m = sprintf('%s\n%s  %s',m,get(hTextClock,'String'),s);
    else
      m = sprintf('%s\n%s  %s',m,get(hTextStateTime,'String'),s);
    end
    if ishandle(hTextMessages)
      set(hTextMessages,'String',m);
      if strncmp(state,'r',1)
        set(hTextMessages,'UserData',m);
      end
    end
    if strcmp(state,'recording') || (strncmp(state,'r',1) && exist([fActive,'messages.mat'],'file'))
      save([fActive,'messages.mat'],'m');
    end
  end

  % saves data to files
  function setSave(full)
    if full || zRawData > 0
      rawdata = dRawData;
      if ~full
        for n = 1:iChannels
          rawdata{n} = rawdata{n}(1:zRawData);
        end
      end
      s = sprintf('rawdata-%03d.mat',getSaveFileNum('rawdata'));
      save([fActive,s],'rawdata');
      setMessage(['Saved raw multichannel data to ',s]);
    end
    if full || zNerves > 0
      for n = 1:iNerves
        nerve = dNerves{n};
        if ~full
          nerve = nerve(1:zNerves);
        end
        s = sprintf('nerve%1d-%03d.mat',n,getSaveFileNum(['nerve',int2str(n)]));
        save([fActive,s],'nerve');
        setMessage(['Saved nerve data to ',s]);
      end
    end
    if iNerves == 2
      nrv2 = tdtRP2;
      nrv2.name = 'Nrv2';
      tdtData = {tdtRX5,tdtRP2,nrv2};
    else
      tdtData = {tdtRX5,tdtRP2};
    end
    if full
      s = 'data.mat';
    else
      s = sprintf('data-%03d.mat',getSaveFileNum('data'));
      if exist([fActive,'data.mat'],'file')
        delete([fActive,'data.mat']);
      end
      for n = 1:iChannels
        dDecData{n} = dDecData{n}(1:zDecData);
      end
      for n = 1:iNerves
        dDecNerves{n} = dDecNerves{n}(1:zDecData);
      end
      % TODO (setSave) append to previous decimated.mat to allow multiple recording sessions
      decimated = struct;
      decimated.dataDecimation = iDecData;
      decimated.nerveDecimation = iDecNerves;
      decimated.data = dDecData;
      decimated.nerves = dDecNerves;
      save([fActive,'decimated.mat'],'decimated');
    end
    save([fActive,s],'tdtData');
    setMessage(['Saved information to ',s]);
    set(hTextFile,'String',int2str(getSaveFileNum('rawdata')));
  end

  % sets spikes in sort graph and last spike for each channel
  function setSpikes(tm,t)
    set(hTextThresh,'String','');
    if get(hGraphThresh,'UIContextMenu') == hMenuThresh
      % TODO (setSpikes) fill tSpikes while editting thresholds (or acquiring), clear when threshold changes
      set(hGraphSpikesLine,'XData',[],'YData',[]);
      d = get(hGraphScope,'UserData');
      if ~isstruct(d) || d.blockUpdate == 0
        return
      end
      for c = 1:iChannels
        if strncmp(getState(),'v',1)
          data = getRawData(c,d.lastFile,max(1,d.lastBlock-d.blockUpdate+1):d.lastBlock);
        else
          % TODO (setSpikes) get proper data while acquiring
          data = [];
        end
        times = getSpikeTimes(data,ySpikes(c),62);
        if c == iChannel
          x = get(hGraphSort,'UserData');
          for spike = getSpikeSnips(data,times,62)
            set(hGraphSortLine(x),'XData',1:length(spike),'YData',spike);
            x = x+1;
            if x > 30
              x = 1;
            end
          end
          set(hGraphSort,'UserData',x);
        else
          spike = [];
          if ~isempty(times)
            spike = data(times(end)-31+(1:62));
          end
        end
        if ~isempty(spike)
          set(hGraphSpikesLine(c),'XData',1:length(spike),'YData',spike);
        end
      end
    elseif ~isempty(dSpikes)
      set([hGraphSortLine,hGraphSpikesLine],'XData',[],'YData',[]);
      for c = 1:iChannels
        r = getRange(tSpikes{c},tm,t);
        if ~isempty(r)
          set(hGraphSpikesLine(c),'XData',1:length(dSpikes{c}(:,r(end))),'YData',dSpikes{c}(:,r(end)));
          if c == iChannel
            if isempty(sSpikes{c})
              setSpikeSorts(c);
            end
            % TODO (setSpikes) display lightgray for sSpikes = 0, iChannel color for sSpikes = 1, iChannel+x for sSpikes = 1+x
            set(hTextThresh,'String',sprintf('%.2f spikes/second (%.0f spikes)',length(r)/(t-tm),length(r)));
            if length(r) > 30
              r = r(end-29:end);
            end
            x = 1;
            for s = r
              set(hGraphSortLine(x),'XData',1:length(dSpikes{c}(:,s)),'YData',dSpikes{c}(:,s));
              x = x+1;
            end
            if x > 30
              x = 1;
            end
            set(hGraphSort,'UserData',x);
          end
        end
      end
    end
  end

  % sets spike shapes for given channel
  function setSpikeShapes(c)
    l = length(tSpikes{c});
    hSpikes{c} = struct;
    hSpikes{c}.time = zeros(1,l);
    hSpikes{c}.height = zeros(1,l);
    for n = 1:l
      [hSpikes{c}.time(n),hSpikes{c}.height(n)] = getSpikeShape(dSpikes{c}(:,n));
    end
  end

  % sets spike sorts for given channel
  function setSpikeSorts(c)
    if isempty(hSpikes{c})
      setSpikeShapes(c);
    end
    sSpikes{c} = ones(1,length(tSpikes{c}));
    ahs = abs(hSpikes{c}.height);
    if get(hCheckSortEasy,'Value') == get(hCheckSortEasy,'Max')
      sSpikes{c}(hSpikes{c}.time < 3 | hSpikes{c}.time > 15 | ahs < 1.2*ySpikes(c) | ahs > 10000) = 0;
    end
    ind = find(sSpikes{c});
    if get(hCheckSortLow,'Value') == get(hCheckSortLow,'Max') && ~isempty(ind)
      ind = find(hSpikes{c}.time > mean(hSpikes{c}.time(ind))+2*std(hSpikes{c}.time(ind)));
      if ~isempty(ind)
        sSpikes{c}(ahs < mean(ahs(ind))+std(ahs(ind))) = 0;
      end
    end
    ind = find(sSpikes{c});
    if get(hCheckSortOutliers,'Value') == get(hCheckSortOutliers,'Max') && ~isempty(ind)
      hmid = mean(ahs(ind));
      hvar = 3*std(ahs(ind));
      tmid = mean(hSpikes{c}.time(ind));
      tvar = 3*std(hSpikes{c}.time(ind));
      sSpikes{c}(ahs < hmid-hvar | ahs > hmid+hvar | hSpikes{c}.time < tmid-tvar | hSpikes{c}.time > tmid+tvar) = 0;
    end
  end

  % sets up TDT interface
  function b = setTDTDefaults(m)
    try
      tdtRX5 = struct;
      tdtRP2 = struct;
      tdtZBus       = tdtGetDevice(hWindow,'zbus.x','ZBUS','GB',1);
      tdtRX5.device = tdtGetDevice(hWindow,'rpco.x','RX5','GB',1);
      tdtRP2.device = tdtGetDevice(hWindow,'rpco.x','RP2','GB',1);
      b = tdtInitializeData(get(hTimer,'UserData'));
    catch
      b = false;
      lasterror
    end
    if b
      cbSetTimeScales(hMenuTimeScale,[]);
    else
      tdtZBus = [];
      tdtRX5  = [];
      tdtRP2  = [];
      if m
        setMessage('Unable to connect to TDT ZBus, RX5, and RP2');
      end
    end
  end

  % returns list of data buffers on a tdt device
  function buffers = tdtGetBuffers(h)
    z = invoke(h,'GetNumOf','ParTag');
    tags = cell(1,z);
    for n = 1:z
      tags{n} = invoke(h,'GetNameOf','ParTag',n);
    end
    tags = unique(tags);
    [par,] = strtok(tags,'~');
    tags = unique(par);
    z = length(tags);
    channels = ones(1,z);
    sizes = ones(1,z);
    data = zeros(1,z);
    for n = 1:z
      temp = tags{n};
      m = length(strmatch(temp,par,'exact'));
      if m > 1
        channels(n) = m;
        if strncmp('d',temp,1)
          sizes(n) = invoke(h,'GetTagSize',temp);
          data(n) = 1;
        end
      elseif invoke(h,'GetTagType',temp) == 68
        sizes(n) = invoke(h,'GetTagSize',temp);
        data(n) = 1;
      end
    end
    nbuffer = 0;
    buffers = cell(1,sum(data));
    for n = find(data)
      nbuffer = nbuffer + 1;
      buffers{nbuffer}.name = tags{n}(2:5);
      buffers{nbuffer}.device = h;
      buffers{nbuffer}.channels = channels(n);
      buffers{nbuffer}.samplingRate = 2*invoke(h,'GetSFreq');
      buffers{nbuffer}.bufferSize = 2*sizes(n);
      buffers{nbuffer}.blockSize = 1;
      buffers{nbuffer}.dataType = 'I16';
    end
    for n = 1:nbuffer
      m = 1;
      f = strfind(tags,buffers{n}.name);
      while m < length(f)
        if ~isempty(f{m})
          switch tags{m}(1)
            case 'b'
              buffers{n}.bufferSize = 2*invoke(h,'GetTagVal',tags{m});
            case 'd'
              buffers{n}.data = tags{m};
            case 's'
              buffers{n}.sync = tags{m};
          end
          f(m) = [];
          tags(m) = [];
        else
          m = m + 1;
        end
      end
    end
  end

  % returns data from a tdt device
  function d = tdtGetData(tdt,tag,offset)
    d = invoke(tdt.device,'ReadTagVEX',tag,offset,tdt.blockSize,tdt.dataType,tdt.dataType,1);
  end

  % returns handle to a tdt device
  function h = tdtGetDevice(figure,activex,device,interface,number)
    h = actxcontrol(activex,[0,0,0,0],figure);
    if strcmp(device,'ZBUS')
      invoke(h,'ConnectZBUS',interface);
    else
      invoke(h,['Connect',device],interface,number);
    end
  end

  % returns status of a tdt device
  function s = tdtGetStatus(h)
    s = invoke(h,'GetStatus');
    if ~bitget(s,1)
      s = 'disconnected';
    elseif ~bitget(s,2) || ~bitget(s,3)
      s = 'stopped';
    else
      s = 'running';
    end
  end

  % returns number of words of data available in data buffer on tdt device
  function r = tdtGetSync(tdt)
    r = invoke(tdt.device,'GetTagVal',tdt.sync);
  end

  % returns processor usage of tdt devices
  function u = tdtGetUsage()
    u = [0;0;0];
    if bitget(invoke(tdtRX5.device,'GetStatus'),3)
      u(1) = invoke(tdtRX5.device,'GetTagVal','zCycUse');
      u(2) = invoke(tdtRX5.device,'GetTagVal','zCycUseAux');
    end
    if bitget(invoke(tdtRP2.device,'GetStatus'),3)
      u(3) = invoke(tdtRP2.device,'GetTagVal','zCycUse');
    end
  end

  % sets up tdt devices
  function b = tdtInitializeData(period)
    set([hGraphNerveLine,hGraphDataLine,hGraphThreshLine,hGraphSpikesLine],'XData',[],'YData',[]);
    set(hGraphScopeLine,'XData',[0,0]);
    set([hTextRX5Data,hTextRP2Data],'String','');
    b = tdtIsReady();
    if b
      tdtSetRun();
      b = bitget(invoke(tdtRX5.device,'GetStatus'),3) && bitget(invoke(tdtRP2.device,'GetStatus'),3);
      if b
        set(get(hToolbar,'Children'),'Enable','off');
        pause(1);
        set(get(hToolbar,'Children'),'Enable','on');
        tdtData = [tdtGetBuffers(tdtRX5.device),tdtGetBuffers(tdtRP2.device)];
        for n = 1:length(tdtData)
          d = invoke(tdtData{n}.device,'GetTagVal',tdtData{n}.sync);
          if d ~= 0
            tdtData{n}.samplingRate = tdtData{n}.samplingRate/round(invoke(tdtData{n}.device,'GetTagVal','zTime')/d);
          end
          tdtData{n}.blockSize = getMaxFactor(tdtData{n}.bufferSize,tdtData{n}.samplingRate*period);
          switch tdtData{n}.name
            case 'Wave'
              tdtRX5 = tdtData{n};
              set(hTextRX5Data,'String',sprintf('%d channels @ %d Hz',tdtRX5.channels,round(tdtRX5.samplingRate)));
            case 'Nrv1'
              tdtRP2 = tdtData{n};
              set(hTextRP2Data,'String',sprintf('1 channel @ %d Hz',round(tdtRP2.samplingRate)));
            case 'Nrv2'
              set(hTextRP2Data,'String',sprintf('2 channels @ %d Hz',round(tdtRP2.samplingRate)));
          end
        end
        tdtRX5.blockSize = tdtRP2.blockSize*tdtRX5.bufferSize/tdtRP2.bufferSize;
        setMessage(sprintf('Data initialized for %d channels and %d nerves',iChannels,iNerves));
      end
      tdtSetIdle();
    end
  end

  % returns whether or not tdt devices are ready
  function b = tdtIsReady()
    b = false;
    if ~isempty(tdtZBus) && isstruct(tdtRX5) && isstruct(tdtRP2)
      try
        b = bitget(invoke(tdtRX5.device,'GetStatus'),1) && bitget(invoke(tdtRP2.device,'GetStatus'),1);
      catch
        lasterror
        tdtZBus = [];
      end
    end
  end

  % set frequencies of band pass filter on tdt rx5 device
  function tdtSetFilterFreqs()
    hp = get(hMenuRX5FilterHighPass,'UserData');
    hp = hp(get(hMenuRX5FilterHighPass,'Value'));
    lp = get(hMenuRX5FilterLowPass,'UserData');
    lp = lp(get(hMenuRX5FilterLowPass,'Value'));
    invoke(tdtRX5.device,'SetTagVal','HPFreq',hp);
    invoke(tdtRX5.device,'SetTagVal','LPFreq',lp);
  end

  % sets tdt devices to idle state
  function tdtSetIdle()
    invoke(tdtZBus,'zBusTrigA',0,2,5);
    invoke(tdtRX5.device,'Halt');
    invoke(tdtRX5.device,'ClearCOF');
    invoke(tdtRP2.device,'Halt');
    invoke(tdtRP2.device,'ClearCOF');
  end

  % sets tdt devices to running state
  function tdtSetRun()
    tdtSetIdle();
    invoke(tdtRX5.device,'LoadCOF',['Multichannel',int2str(iChannels),'.rcx']);
    invoke(tdtRX5.device,'Run');
    tdtSetFilterFreqs();
    invoke(tdtRP2.device,'LoadCOF',['Nerves',int2str(iNerves),'.rcx']);
    invoke(tdtRP2.device,'Run');
    invoke(tdtZBus,'zBusTrigA',0,1,5);
  end

end
