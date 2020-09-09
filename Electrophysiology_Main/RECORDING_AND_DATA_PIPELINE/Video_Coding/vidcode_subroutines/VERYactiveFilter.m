function [datback] = VERYactiveFilter (mvmt, pts, tX, dat)
% Now redo the same thing but for a VERY HIGH movement threshold. this
% should override anything else - artifact can mess up the EMG, LFP, etc,
% but total video movement is always clear when reviewed. This is the best
% code based approximation of that rule.
% find crossings of AWTHRESH line - these seem to correspond with active wake:

Vawidx(:,1) = mvmt(:,2);
Vawidx(:,2) = mvmt(:,1) >= 1.35;


vadiff = diff(Vawidx(:,2));
vawons = find(vadiff == 1) + 1;
vaoffs = find(vadiff == -1) + 1;

if size(vawons,1) ~= size(vaoffs,1);
    
    if Vawidx(1,2) == 1;
        
        vawons = [1; vawons]; % if it starts with AW
        
    elseif Vawidx(end,2) == 1;
        
        vaoffs = [vaoffs; numel(Vawidx(:,2))]; % if it ends with AW
        
    end
end

for ee = 1:size(vawons,1);
    VAWIDX(ee,1) = vawons(ee);
    VAWIDX(ee,2) = vaoffs(ee);
    
    VAWTIMES(ee,1) = Vawidx(vawons(ee));
    VAWTIMES(ee,2) = Vawidx(vaoffs(ee));
    
end


if ~isempty(vawons);
    for ee = 1:size(VAWTIMES,1);
        pts(tX >= VAWTIMES(ee,1) & tX <= VAWTIMES(ee,2),2 ) = 4 ; % Be consistent. 4 is the standard code used for active waking in our software
    end
end


% High EMG can also be a good sign of artifact in the LFP trace - consider
% a case when animals are rough housing and generate false high
% delta/theta. this should appear in movement (video) as well as in the
% EMG. 

hEMG = (dat(:,3)>1.5);
temp = pts(:,2)<3;
pts(hEMG&temp,2) = 4;


% Finally, clean up points that are inappropriately coded as active. look
% for the intersection of wake coded points and high delta (above some
% threshold. at this point [11/15] we have a lot of artifact in the video
% trace, so it's not a perfectly reliable indicator of being awake)


awtmp       = pts(:,2)>3;
ttemp       = dat(:,2)>1.25;
dtemp       = dat(:,1)>1.5;
emgtemp     = dat(:,3)<0.8;
collusion   = awtmp & dtemp & ttemp & emgtemp;

pts(collusion,2) = 2;

% Test this function: look for instances where delta and theta are
% intermediate, and there IS movement. Call this quiet wake.

stmp    = pts(:,2)<3;
ttemp   = dat(:,2)<1.5;
dtemp   = dat(:,1)<2;

qwtemp = stmp & ttemp & dtemp;
qd = diff(qwtemp);
qON = find(qd == 1);
qOFF = find(qd == -1);

if qwtemp(1,1) == 1;
    qON = [1; qON];
end

if qwtemp(end,1) == 1;
   qOFF = [qOFF; size(pts,1)]; 
end

qIDX = [qON qOFF];

qIDX(qIDX(:,2) - qIDX(:,1)<100,:) = [];
timeQ = tX(qIDX);

if ~isempty(timeQ);
    for ee = 1:size(timeQ,1);
        mcheck(ee) = mean ( mvmt(  mvmt(:,2)>=timeQ(ee,1) & mvmt(:,2)<=timeQ(ee,2), 1 ) );
    end
    
    
    mcheck = (mcheck > 0.6)'; % set some threshold for movement that must occur during ambiguous quiet waking.
    
    qIDX = qIDX(mcheck,:);
    
    for ee = 1:size(qIDX,1);
        
        pts(qIDX(ee,1):qIDX(ee,2),2) = 5;
        
    end
end

datback = pts;





