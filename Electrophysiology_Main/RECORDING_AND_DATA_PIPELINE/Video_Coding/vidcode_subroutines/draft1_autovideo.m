%
%
%       NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5
%
%

clearvars -except LFPinfo outdata sd_save EMGdata
close all;

addpath('/Users/khengen/Desktop/VIDEOSTART/');

if ~exist('LFPinfo','var');
    disp('loading lfp info');
load('/Users/khengen/Desktop/VIDEOSTART/LFPinfo.mat');
end
if ~exist('sd_save','var');
    disp('loading SW data');
load('/Users/khengen/Desktop/VIDEOSTART/KH57_SD_states.mat');
end
if ~exist('outdata','var');
    disp('loading movement information');
    load('/Users/khengen/Desktop/VIDEOSTART/KH57_movement_v2.mat');
end
if ~exist('EMGdata','var');
    [eFile, eDir] = uigetfile; % load emg file
    load([eDir eFile]);
end



BLOCK = 7;

% LFP start time is NOT GMT corrected - do this.
if ~isfield(LFPinfo,'GMTcor');
    for ee = 1:size(LFPinfo,2);
        LFPinfo(ee).startTime    =  LFPinfo(ee).startTime - 3600*4;
        LFPinfo(ee).GMTcor       = 1;
    end
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Eventually, start the looping here:





P = LFPinfo(BLOCK).spectrogram.P;
F = LFPinfo(BLOCK).spectrogram.F;

pX = linspace(LFPinfo(BLOCK).startTime,  LFPinfo(BLOCK).startTime + LFPinfo(BLOCK).duration, size(P,2) );

deltapower  = sum(10*log10(abs(P(F<=4,:))));
deltapower  = deltapower/abs(max(deltapower));
deltapower  = 10*(deltapower+abs(min(deltapower))); % make positive and scale for the plot
deltapower  = deltapower.^1.5;
deltapower  = smooth(deltapower,100);

thetapower  = sum(10*log10(abs(P(F>=6 & F<=8,:))));
thetapower  = thetapower/abs(max(thetapower));
thetapower  = 10*(thetapower+abs(min(thetapower)));
thetapower  = smooth(thetapower,200);
thetapower  = thetapower.^2;





emgX = linspace(LFPinfo(BLOCK).startTime,  LFPinfo(BLOCK).startTime + LFPinfo(BLOCK).duration, size(EMGdata.f_EMG_L,2) );       
%emgX = %emgX * (size(P,2)/size(EMGdata.f_EMG_L,2));
emgY = 10.*(EMGdata.f_EMG_L/max(EMGdata.f_EMG_L));

semgY = smooth(EMGdata.f_EMG_L,2000);
semgY = 10.*(semgY/max(semgY));
semgY = semgY/mean(semgY);


% take a moving average to smooth/resample the EMG
bb = round(linspace(1, size(semgY,1), size(pX,2)+1  ));

for ee = 1:size(bb,2)-1;
    rs_emg(ee) = mean( semgY( bb(ee):bb(ee+1) ));
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% figure out the state coded plot here:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
statetimes = sd_save(:,1:2);
statetimes = sortrows(statetimes,2);

testD   = diff(statetimes(:,1));
killem  = find(testD == 0);
killem  = killem+1;
statetimes(killem,:) = [];

% get rid of 3's: general wake
threes = [];
threes = find(statetimes(:,1) == 3);

if ~isempty(threes);
    
    if threes(end) == size(statetimes,1);
        statetimes(threes(1:end-1),1) = statetimes(threes(1:end-1)+1,1);
        statetimes(end,1) = 4;
    else
        
        statetimes(threes,1) = statetimes(threes+1,1);
    end
    
    % Clean up multiple button presses for a single bout
    testD   = diff(statetimes(:,1));
    killem  = find(testD == 0);
    killem  = killem+1;
    statetimes(killem,:) = [];
end

s1 = find(statetimes(:,2) >= LFPinfo(BLOCK).startTime,1); % this is the first state transition greater than the LFP start time
s2 = find(statetimes(:,2) >= LFPinfo(BLOCK).startTime + LFPinfo(BLOCK).duration,1);

tState = statetimes(s1:s2,:);

if statetimes(s1,2) < LFPinfo(BLOCK).startTime;
    tState(1,2) = LFPinfo(BLOCK).startTime;
end

if statetimes(s2,2) > LFPinfo(BLOCK).startTime + LFPinfo(BLOCK).duration;
   tState(end,2) =  LFPinfo(BLOCK).startTime + LFPinfo(BLOCK).duration;
end

awAREA = find(tState(1:end-1,1) == 4);
qwAREA = find(tState(1:end-1,1) == 5);
rmAREA = find(tState(1:end-1,1) == 1);
nrAREA = find(tState(1:end-1,1) == 2);



d_power = deltapower>1.5;
t_power = thetapower>0.5;



% Process the video analysis data for plotting:
mvmt = outdata.movement;
mvmt(mvmt(:,3)<pX(1),:) = [];
mvmt(mvmt(:,3)>pX(end),:) = [];
mvmt(:,2) = smooth(mvmt(:,2),60);
mvmt(:,2) = mvmt(:,2)/mean( mvmt( mvmt(:,2)~=0 ,2)   );


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - AUTOMATE THE EXTRACTION PROCESS - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % put the relevant data into a single matrix for easy passing between
% functions
sm_delta = smooth(deltapower,20);
seq = [sm_delta thetapower rs_emg'];


% set the minimum epoch time - this is the threshold for analysis of
% potential errors based on duration
T_threshWAKE = 60;
T_threshSLEEP = 30; % also includes REM

% find waking states
WAKEverbose = 0;
w_dat   = auto_wake(seq, pX, T_threshWAKE, mvmt,WAKEverbose); 

% find NREM and REM states
s_dat = auto_sleepStates(seq, pX, T_threshSLEEP, w_dat);
% START HERE - make sleep states ignore the waking states that are already
% encoded
% figure(77)
% plot(s_dat(:,1),s_dat(:,2) == 1)
% set(gca,'ylim',[0 4]);

[datback]   = VERYactiveFilter (mvmt, s_dat, pX);

[datback1]  = QWREMfilter (mvmt, seq, datback, pX) ;

minT = 60;
[datback2]  = flicker(datback1,minT,pX);

% - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% edit out REM that occurs during inappropriate blocks:

[datback3]  = REMeditor(datback2, seq);
% figure(78)
% plot(datback3(:,1),datback3(:,2) == 1)
% axis([datback3(1,1) datback3(end,1) 0 4]);




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




% COMPARE THIS VISUALLY:
stse = [];
stse(:,1) = pX;
stse(:,2) = datback3(:,2);

tD      = diff(stse(:,2));
tD(end) = 1;
killem  = find(tD == 0);
killem  = killem+1;
stse(killem,:) = [];

%stse = [stse; [pX(end) 0] ];

AUTrmAREA       = find(stse(1:end-1,2) == 1);
AUTnrAREA       = find(stse(1:end-1,2) == 2);
AUTawakeAREA    = find(stse(1:end-1,2) == 4);
AUTqwakeAREA    = find(stse(1:end-1,2) == 5);




% - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - 
%  - - - - - - - - - - - DO THE PLOTTING - - - - - - - - - - - - -
% - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - 

% Plot the spectrogram:
%close all;

figure(886);

imagesc(pX,( F ),10*log10(abs(P)));
set(gca,'ydir','normal')
hold on
plot(pX, thetapower,'b','linewidth',2);
plot(pX, sm_delta,'r','linewidth',2);
plot(mvmt(:,3),mvmt(:,2),'m');
%plot(pX, d_power,'r');
%plot(pX, t_power,'b');

% Plot the EMG trace:

plot(pX, rs_emg);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot the manually coded video:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for ee = 1:size(awAREA,1)
    
    rectangle( 'Position', [tState(awAREA(ee),2), 15, tState(awAREA(ee)+1,2) - tState(awAREA(ee),2) , 1],...
        'facecolor','g','linestyle','none'), hold on
    
end
for ee = 1:size(qwAREA,1)
    
    rectangle( 'Position', [tState(qwAREA(ee),2), 15, tState(qwAREA(ee)+1,2) - tState(qwAREA(ee),2) , 1],...
        'facecolor','y','linestyle','none'), hold on
    
end
for ee = 1:size(rmAREA,1)
    
    rectangle( 'Position', [tState(rmAREA(ee),2), 15, tState(rmAREA(ee)+1,2) - tState(rmAREA(ee),2) , 1],...
        'facecolor',[0.5 0.5 0.5],'linestyle','none'), hold on
    
end
for ee = 1:size(nrAREA,1)
    
    rectangle( 'Position', [tState(nrAREA(ee),2), 15, tState(nrAREA(ee)+1,2) - tState(nrAREA(ee),2) , 1],...
        'facecolor',[0.2 0.2 0.3],'linestyle','none'), hold on
    
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot the algorithmically coded video:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

for ee = 1:size(AUTawakeAREA,1)
    
   
        rectangle( 'Position', [stse(AUTawakeAREA(ee),1), 16, stse(AUTawakeAREA(ee)+1,1) - stse(AUTawakeAREA(ee),1) , 1],...
            'facecolor','g','linestyle','none'), hold on
  
end

for ee = 1:size(AUTqwakeAREA,1)
    
    rectangle( 'Position', [stse(AUTqwakeAREA(ee),1), 16, stse(AUTqwakeAREA(ee)+1,1) - stse(AUTqwakeAREA(ee),1) , 1],...
        'facecolor','y','linestyle','none'), hold on
    
end

for ee = 1:size(AUTrmAREA,1)
    
    rectangle( 'Position', [stse(AUTrmAREA(ee),1), 16, stse(AUTrmAREA(ee)+1,1) - stse(AUTrmAREA(ee),1) , 1],...
        'facecolor',[0.5 0.5 0.5],'linestyle','none'), hold on
    
end
for ee = 1:size(AUTnrAREA,1)
    
    rectangle( 'Position', [stse(AUTnrAREA(ee),1), 16, stse(AUTnrAREA(ee)+1,1) - stse(AUTnrAREA(ee),1) , 1],...
        'facecolor',[0.2 0.2 0.3],'linestyle','none'), hold on
    
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


set(gca,'ylim',[0.25 17]);
set(gca,'TickLength',[0 0.025]);







% DEAL WITH BOOKENDS - EXTEND THE FIRST AND LAST MINUTE?
