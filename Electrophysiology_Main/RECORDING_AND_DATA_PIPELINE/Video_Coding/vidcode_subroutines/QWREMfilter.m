function [dataout] = QWREMfilter (~, seq, datain, ~)
%
% [DATAOUT] = QREMFILTER(MVMT, SEQ, DATAIN, PX)
%
% Evaluate all instances of quiet waking and check to see if they might
% actually be better explained by REM. This is established by checking the
% variance of the EMG trace and the prior behavioral state. If the prior
% state was a sleep state and there is VERY low EMG variance, the epoch of
% quiet waking will be over-written as REM sleep.
%
% Inputs:
% SEQ: The first column is the integral of the delta band power, the second
% column is the integral of the theta band power, and the third column is
% the EMG integral.
% DATAIN: n points x 2 matrix. The first column is the sample point index
% and the second column is the animal behavioral state code.
% PX: time stamps that correspond to the sample index in DATAIN(:,1).
% MVMT: Movement data, sampled at a different rate than SEQ, but with
% timestamps that can be aligned using PX. This and PX are included for
% future use of MVMT to disambiguate REM and quiet waking should the need
% arise.
%
% Output:
%   DATAOUT is a n samples x 2 matrix. The first column contains the time
%   indices of sampling. The second column contains the categorical output
%   of Active Wake (4), Quiet Wake (5), REM (1), or NREM (2).
%
%   11/2015 KBH


statetimes = datain;

% find state transitions:
testD   = diff(statetimes(:,2));
killem  = find(testD == 0);
killem  = killem+1;
statetimes (killem,:) = [];

q = find(statetimes(:,2) == 5);

if ~isempty(q); % make sure there are any quiet waking bouts in this block
    
    QWon = statetimes(q,1);
    
    
    if QWon(end) == statetimes(end,1);
        
        QWoff= statetimes(q(1:end-1)+1 ,1);
        QWoff = [QWoff; datain(end,1)];
    else
        QWoff= statetimes(q+1 ,1);
    end
    
    
    
    if size(QWon,1) ~= size(QWoff,1);
        disp('Deal with mismatched index sizes in QWREMfilter.m')
        keyboard
    end
    
    
    % find quiet wake segments with really, really low variance in the EMG
    % trace. If these periods were preceeded immediately by a sleep state,
    % they're probably REM and not QW.
    
    for ee = 1:size(QWon);
        
        qvar(ee) = var(seq(QWon(ee,1):QWoff(ee,1),3));
        
        if qvar(ee)<0.01;
            
            if datain(QWon(ee)-1,2)<3;
                
                datain(QWon(ee):QWoff(ee)-1,2) = 1;
                
            end
            
        end
        
    end
    
    % Find points that are called QW that also have low delta and high theta.
    % Re evaluate these as REM states instead.
    qwtemp  = datain(:,2) == 5;
    thta    = seq(:,2)>1.5;
    delta   = seq(:,1)<1.5;
    emg     = seq(:,2)<0.5;
    tempREM = qwtemp & thta & delta & emg;
    datain(tempREM,2) = 1;
    
end

dataout = datain;


