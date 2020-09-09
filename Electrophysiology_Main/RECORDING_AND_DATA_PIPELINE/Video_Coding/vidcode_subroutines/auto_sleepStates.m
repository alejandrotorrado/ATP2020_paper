function [sleep] = auto_sleepStates(dat, tpts, minT, wake)
%
% [SLEEP] = AUTO_SLEEPSTATES(DAT, TPTS, MINT, WAKE)
%
% Return a categorical assessment of REM and NREM states for each sample
% point passed to the function in DAT. Data are a set of polysomnographic 
% measures including the integrated delta and theta power of the LFP, as 
% well as the integrated EMG signal recorded simultaneously. DATA should be
% organized as n sample points x 3 columns.
%
% The first column is the integral of the delta band power, the second
% column is the integral of the theta band power, and the third column is
% the EMG integral. TPTS is the 1 x nsample points vector containing the 
% time indices of each data points in DAT. MINT is the minimum time that 
% an epoch will be considered valid. WAKE is the previously extracted
% matrix of waking times. These are automatically excluded from the
% analysis of sleep states (N.B. REM can be difficult to distinguish from
% waking without an EMG signal. EMG is used in AUTO_WAKE to determine when
% an animal is awake and asleep).
%
%   NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5
%
% Output:
%   SLEEP is a n samples x 2 matrix. The first column contains the time
%   indices of sampling. The second column contains the categorical output
%   of wake (NaN), REM (1), or NREM (2). 
%
%   10/2015 KBH

pts(:,1)    = 1:size(tpts,2);
pts(:,2)    = wake(:,2);


dataidx(:,1) = 1:length(tpts);
dataidx(:,2) = (dat(:,1) - dat(:,2))>0.2 & wake(:,2)~=1 ; % find the NREM samples
dataidx(:,2) = dataidx(:,2)*2; % NREM is a "2" in our body of code. Keep it consistent.
dataidx(wake(:,2)~=0,2) = wake(wake(:,2)~=0,2); 
dataidx(dataidx(:,2) == 0,2) = 1; % The remaining bins are, by definition, REM (assuming everything is working properly). REM is a "1"

%dataidx(isnan(dataidx(:,2)), 2) = 3; % general wake states should be "3"

% Look for inappropriate state transitions and correct them:
sts(:,1) = tpts;
sts(:,2) = dataidx(:,2);

tD      = diff(sts(:,2));
killem  = find(tD == 0);
killem  = killem+1;
sts(killem,:) = [];

tmp(:,1) = sts(1:end-1,2);
tmp(:,2) = sts(2:end,2);


% Look for periods of high delta that are persistant and then make sure
% these are tagged as NREM:

DTHR    = 1.85;
TTHR    = 1.5;
nrm     = dat(:,1)>DTHR; % find delta datapoints above threshold
theta   = dat(:,2)>TTHR;
nrem    = nrm & theta;

dataidx(nrem,2) = 2;

% dnrem   = diff(nrem); % find onset and offset of high delta periods
% nremON  = find(dnrem == 1); % indices of the onsets
% nremOFF = find(dnrem == -1); % indices off the offsets


% 
% if dat(1,1)>DTHR;
%     nremON = [1; nremON];
% end
% 
% if dat (end,1)>DTHR;
%     nremOFF = [nremOFF; numel(dat(:,1))];
% end
% 
% nrbout  = [nremON nremOFF];
% nrdur   = tpts(nrbout);
% nrdur   = nrdur(:,2) - nrdur(:,1);
% 
% keyboard
% 
% for ee = size(nremON,1);
%     dataidx(nremON(ee):nremOFF(ee),2) = 2;
% end

% 
% figure(7)
% plot(tpts,dat(:,1)), hold on
% plot(tpts,dataidx(:,2) == 2)
% 
% plot(tpts,dat(:,1))
% for ee = 1:size(nremOFF,1);
% 
%     e(ee) = line([tpts(nremON(ee,1)) tpts(nremON(ee,1))], [0 4]);
%     set(e(ee),'color','g');
%     f(ee) = line([tpts(nremOFF(ee,1)) tpts(nremOFF(ee,1))], [0 4]);
%     set(f(ee),'color','r');
% 
% end
% 
% axis([tpts(1) tpts(end) 0 4]);
% hold off


% flag wake to REM transitions:

a       = dataidx(1:end-1,2);
a(:,2)  = dataidx(2:end,2);
da      = [0; diff(a(:,1))];
w2r     = find(a(:,1)>3 & a(:,2) == 1);

dataidx = [dataidx; dataidx(end,1)+1 NaN];
if ~isempty(w2r);
    disp('Correcting REM immediately following wake.') 
    for ee = 1:length(w2r);
        disp([num2str(ee) ' of ' num2str(length(w2r)) '.']);
        
        kstart  = w2r(ee)+1 ; % find the error REM "on" sample
        kend    = kstart + find(da(kstart:end)~=0,1); % find the error REM "off" sample
        
        
        if dataidx(kstart-1,2) < 3 | dataidx(kend+1,2) < 3;
            % if the REM is flanked on either side by NREM, convert the
            % chunk to NREM
            dataidx(kstart:kend,2) = 2;
            
        else
            % if it's flanked by waking on either side, convert chunk to QW
            dataidx(kstart:kend,2) = 4; % set the error sample points to " quiet wake"

        end
        
    end

end
    
% 
% statecode = 1; % rem
% %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% datback = flicker(statecode, dataidx,tpts,minT);
% 
% %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% statecode = 2; % nrem
% %  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% datback = flicker(statecode, datback,tpts,minT);

%  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if any(isnan(dataidx(end,:)));
    dataidx(end,:) = [];
end


sleep = dataidx;


% check the data:
%    figure(7), plot(tpts, dataidx(:,2)), hold on
%    plot(tpts, dat(:,1));

