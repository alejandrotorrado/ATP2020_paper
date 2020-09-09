function [events] = GetEventNames (Tank, Block, varargin)
%This function returns all of the events in the specified block.  The
%arguments are as follows:
%Tank is the name of a registered Tank, or the full path with name of an
%unregistered Tank, ie 'DEMOTANK2' or 'C:\mytanks\tank1'
%Block is simply the blockname, ie 'Block-1'
%The final optional argument is for the event type; allowed types are
%'StrobeOn', 'StrobeOff', 'Scalar', 'Stream', 'Snip', 'Mark', and
%'HasData'.  Please see the OpenDeveloper manual for an explanation of the
%different event types, which can be downloaded here: 
%http://www.tdt.com/T2Download/manuals/OpenDeveloper_Manual.pdf

numvarargs = length(varargin);
if numvarargs > 1
    error('Too many arguments');
elseif numvarargs == 0
    type = 'skip';
elseif numvarargs == 1
    type = varargin{1};
end

if strcmp(type,'StrobeOn')
    EventCode = 257;
elseif strcmp(type,'StrobeOff')
    EventCode = 258;
elseif strcmp(type,'Scalar')
    EventCode = 513;
elseif strcmp(type,'Stream')
    EventCode = 33025;
elseif strcmp(type,'Snip')
    EventCode = 33281;
elseif strcmp(type,'Mark')
    EventCode = 34817;
elseif strcmp(type,'HasData')
    EventCode = 32768;
elseif strcmp(type,'skip')
    EventCode = 0;    
else
    error('Event type not understood, allowed types are StrobeOn, StrobeOff, Scalar, Stream, Snip, Mark, and HasData');
end

events = [ ];
TTX = actxcontrol('TTank.X');
if TTX.ConnectServer('Local','Me') == 0 error('Error connecting to server'); end
if TTX.OpenTank(Tank,'R') == 0 error('Error opening tank'); end
if TTX.SelectBlock(Block) == 0 error('Error opening block'); end

codes = TTX.GetEventCodes(EventCode);
num = length(codes);

if isnan(codes) ~= 1
    strArray = java_array('java.lang.String', num);

    for i = 1:num
        strArray(i) = java.lang.String(TTX.CodeToString(codes(i)));
    end
        
    events = cell(strArray);  
    
end
    
if isempty(events)
    clear events;
    fprintf('No events with entered Event Code\n');
end

end



