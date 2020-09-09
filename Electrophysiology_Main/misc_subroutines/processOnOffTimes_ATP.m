function [outrt,outT] = processOnOffTimes_ATP(FR,TLIST,binsz,dayExpStart,dataType,verbose)
%
% outrt = processOnOffTimes(FR,TLIST,binsz)
%
% INPUTS:
% __FR:     should be a column vector containing the firing rate values.
% __TLIST:  should be a list of on/off times, with all on times in the
%           first column and all off times in the second.
% __binsz:  the size (in seconds) of the bin used to calculate the FR
% __dayExpStart: experiment start day. Will adjust on/off times to this.
% __dataType: dataset type ('ER','recov','CPP', etc.). Will determine how to
%             deal with misaligned on/off times.
%
% OUTPUT:
% __outrt:  column vector with firing rate values. The values between the
%           on/off times in TLIST will have been converted to NaNs.
% __outT: list of on/off times
%
% ATP 12/12/2015
%
%
% keyboard;


%% setup
if nargin == 5
    verbose = 1;
end

% depending on dataset, change "check" day
switch dataType
    case {'ER'}
        dayCheck = dayExpStart;
    case {'recov','CPPrecov','CPP2recov','SD_NEW'}
        dayCheck = 6;
end

%% adjust on/off times if needed
% if ontimes are less than the day above, they are not adjusted to BL1
if TLIST(1)/(24*3600) < dayCheck
    if verbose, disp('Adjusting on/off times'); end
    TLIST = TLIST + dayExpStart*24*3600;
else
    if verbose, disp('On/Off times look already adjusted. Make sure!'); end
end

%% apply on/off times to FR array

% case of one pair of on/off times
if size(TLIST,1) == 1 % one on/off time
    tstart  = ceil(TLIST(1,1)/binsz);
    tstop   = floor(TLIST(1,2)/binsz);
    
    if tstart == 0
        tstart = 1;
    end
    if tstop > size(FR,1)
        tstop = size(FR,1);
    end
    if tstop == 0
        tstop = 1;
    end
    
    FR(1:tstart,1)  = NaN;
    FR(tstop:end,1) = NaN;
    
    clear tstart tstop

% case with multiple on/off times
else
    if verbose, fprintf('Multiple on/off times found.\n'); end
    
    for xx = 1:size(TLIST,1)
        tstart  = floor(TLIST(xx,1)/binsz);
        tstop   = ceil(TLIST(xx,2)/binsz);
        
        if tstart == 0
            tstart = 1;
        end
        if tstop > size(FR,1)
            tstop = size(FR,1);
        end
        
        if xx == size(TLIST,1)
            tlast = ceil(TLIST(xx-1,2)/binsz);
            FR(tlast:tstart,1) = NaN;
            FR(tstop:end,1) = NaN;
            clear tlast
        else
            if xx > 1
                tlast = ceil(TLIST(xx-1,2)/binsz);
            else
                tlast = 1;
            end
            tnext = ceil(TLIST(xx+1,1)/binsz);
            if tnext > size(FR,1)
                tnext = size(FR,1);
            end
            FR(tlast:tstart,1) = NaN;
            FR(tstop:tnext,1) = NaN;
            clear tnext tlast
        end
        
        clear tstart tstop 
    end
end

outrt = FR;
outT =  TLIST; % optional output

