function [wakeID] = auto_wake(data, tX, Tmin, mvmt, verbose)
%
% [WAKEID] = AUTO_WAKE(DATA, DELTA, TX)
%
% Return a binary assessment of waking or not waking for each sample point
% of a set of polysomnographic data including the integrated delta and
% theta power of the LFP, as well as the integrated EMG signal recorded
% simultaneously. DATA should be organized as n sample points x 3 columns.
% The first column is the integral of the delta band power, the second
% column is the integral of the theta band power, and the third column is
% the EMG integral. TX is the 1 x nsample points vector containing the time
% indices of each data points in DATA. TMIN is the minimum time that an
% epoch will be considered valid. Periods of wake or not wake that have a
% duration of less than TMIN will be considered for elimination.
%
%       NOTE: Organization of DATA = [deltapower thetapower rs_emg'];
%       NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5
%
% Output:
%   WAKEID is a n samples x 2 matrix. The first column contains the time
%   indices of sampling. The second column contains the binary decision of
%   wake (1) or not wake (0).
%

% create the sample time and state assessment matrix that will hold all
% relevant information for output:


if ~isempty (mvmt) && ~isnan(sum(mvmt(:,2))) ;
    pts(:,1)    = 1:size(tX,2);
    pts(:,2)    = zeros;
    
    % Process the video analysis data:
    % mvmt = vidtrack.movement;
    % mvmt(mvmt(:,3)<tX(1),:) = [];
    % mvmt(mvmt(:,3)>tX(end),:) = [];
    % mvmt(:,2) = smooth(mvmt(:,2),60);
    % mvmt(:,2) = mvmt(:,2)/mean( mvmt( mvmt(:,2)~=0 ,2)   );
    
    AWTHRESH = 0.6;
    
    % find crossings of AWTHRESH line - these seem to correspond with active wake:
    awidx(:,1) = mvmt(:,2);
    awidx(:,2) = mvmt(:,1) >= AWTHRESH;
    
    adiff = diff(awidx(:,2));
    awons = find(adiff == 1) + 1;
    aoffs = find(adiff == -1) + 1;
    
    if size(awons,1) ~= size(aoffs,1);
        
        if awidx(1,2) == 1;
            
            awons = [1; awons]; % if it starts with AW
            
        elseif awidx(end,2) == 1;
            
            aoffs = [aoffs; numel(awidx(:,2))]; % if it ends with AW
            
        end
    end
    
    % get the on and off time of the movement-based AW indexing
    for ee = 1:size(awons,1);
        AWIDX(ee,1) = awons(ee);
        AWIDX(ee,2) = aoffs(ee);
        
        AWTIMES(ee,1) = awidx(awons(ee));
        AWTIMES(ee,2) = awidx(aoffs(ee));
        
    end
    

    % Plot this if you want to check the relative impact of each variable:
    
    if verbose == 1;
        figure (3);
        plot(tX,data(:,1)),hold on
        area(tX,data(:,1)>1.5); % binary decision of NREM sleep?
        plot(mvmt(:,2),mvmt(:,1),'m')
        plot(tX,data(:,3),'--k');
        plot(tX,data(:,2),'b','linewidth',2);
        p = line([tX(1) tX(end)],[AWTHRESH AWTHRESH]);
        set(p,'color','r','linestyle','--');
        xlim([tX(1) tX(end)]);
        legend('Delta','Binary NREM','Movement','EMG','Theta',[num2str(AWTHRESH) ' Thresh'])
        title('All traces during search for waking');
        
        for uu = 1:size(AWIDX,1);
            
            f(uu)    = line([AWTIMES(uu,1) AWTIMES(uu,1)], [0 4 ] );
            set(f(uu),'color','g');
            stf(uu) = line([AWTIMES(uu,2)  AWTIMES(uu,2)], [ 0 4] );
            set(stf(uu),'color','r')
            
        end
        
    end
    
    
    
    % Code active waking based on the AWTHRESH but NOT IF DELTA IS HIGH!
    for ee = 1:size(AWTIMES,1);
        pts(tX >= AWTIMES(ee,1) & tX <= AWTIMES(ee,2),2 )  = 4 ; % & data(:,1)<2.5 Be consistent. 4 is the standard code used for active waking in our software
    end
    
    pts( pts(:,2) == 4 & data(:,1)>1.75 & data(:,2)>0.75, 2 ) = 2;
    
    
    
    % Look for crossings of the QW thresh - this will be used in combination
    % with a few other variables to find waking periods that are not associated
    % with a lot of obvious movement.
    QWTHRESH = 0.5;
    qwidx(:,1) = mvmt(:,2);
    qwidx(:,2) = mvmt(:,1)>=QWTHRESH & mvmt(:,1)<AWTHRESH;
    
    qdiff = diff(qwidx(:,2));
    qwons = find(qdiff == 1) + 1;
    qoffs = find(qdiff == -1) + 1;
    
    if size(qwons,1) ~= size(qoffs,1);
        
        if qwidx(1,2) == 1;
            qwons = [1; qwons];
        end
        
        if qwidx(end,2) == 1;
            qoffs = [qoffs; numel(qwidx(:,2))];
        end
        
    end
    
    for ee = 1:size(qwons,1);
        
        qWIDX(ee,1) = qwons(ee);
        qWIDX(ee,2) = qoffs(ee);
        
        QWTIMES(ee,1) = qwidx(qwons(ee));
        QWTIMES(ee,2) = qwidx(qoffs(ee));
        
    end
    
    qtest(:,1) = tX;
    for ee = 1:size(QWTIMES,1);
        qtest(tX >= QWTIMES(ee,1) & tX <= QWTIMES(ee,2),2) = 1;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    % You'll have to do something to compare the temp value below with an
    % intermediate movement trace, maybe around 0.5, to get the quiet waking.
    % Maybe run temp first, then just assign the remaining time (that's not AW
    % but still general wake) as quiet. You'll have to dig into places where
    % the two measures give opposite conclusions.
    
    
    
    temp    = data(:,3)>data(:,2) & data(:,1)<0.7*mean(data(:,1));
    temp2   = data(:,1)<1.75 & qtest(:,2);
    
    %%%%
    
    if verbose == 1;
        figure (1), hold on
        area(tX,temp2,'facecolor','g');
        axis([tX(1) tX(end) 0 4]);
        title('Inital search for quiet waking')
    end
    
    pts(temp2,2) = 5; % be consistent. Quiet waking is coded as "5" in our software.
    
    
    % reconcile TEMP and the current AW and QW standings. TEMP uses EMG to evaluate wake and sleep
    
    % need to set it so that if there is NO movement coming through the video
    % but there is good EMG data, the EMG will be used in place of blob
    % analysis
    
    % start by assigning EMG identified regions as 3's - convert these to quiet
    % wake eventually once they're confirmed
    pts(pts(:,2)==0 & temp,2) = 3;
    
    
    % now delete 3's that occur in isolation assuming there is otherwise video
    % data used in this block
    if any( unique(pts(:,2)) == 4 | unique(pts(:,2)) == 5 );
        
        % find the onset and offset of "3" states
        idx3    = pts(:,2) == 3;
        d3      = diff(idx3);
        on3     = find(d3 == 1) + 1;
        off3    = find(d3 == -1) + 1;
        
        % Correct for case in which the block ends with a "3" state
        if size(on3,1) > size(off3,1) && idx3(end) == 1;
            off3 = [off3; numel(idx3)];
        end
        
        
        for ee = 1:size(on3,1);
            
            if on3(ee) - 20 < 0;
                prior = 0;
            else
                prior   = mode( pts(on3(ee)-20 : on3(ee),2) );
            end
            
            if off3(ee) + 20 > numel(pts(:,2));
                post = 0;
            else
                post    = mode( pts(off3(ee) : off3(ee)+20,2) );
            end
            
            pvals = [prior post];
            
            % case in which the "3" is surrounded by nothing.
            if post == 0 && prior == 0;
                pts(on3(ee):off3(ee),2) = 0;
            end
            
            % case in which "3" is flanked by a zero on one side and a "4" or
            % "5" on the other side.
            if any(pvals == 0) && post~=prior
                pts(on3(ee):off3(ee),2) = pvals(pvals~=0);
            end
            
            % case in which "3" is flanked on each side by incongruent non zero elements.
            if all(pvals~=0) && prior~=post;
                halfway = round(on3(ee) + (off3(ee)-1 - on3(ee))/2 );
                pts(on3(ee):halfway)    = prior;
                
                pts(halfway+1:off3(ee),2) = post;
                
            end
            
        end
        
        
    end
    
    
    
    % if theta is high and movement low
    combo = []; t = []; td = []; ton = []; toff = []; tms = []; thta = [];
    
    combo(:,1) = tX;
    t       = mvmt(:,1)<0.75;
    td      = diff(t);
    ton     = find(td == 1);
    toff    = find(td == -1);
    if t(1) == 1, ton = [1; ton]; end
    if t(end) == 1; toff = [toff; numel(t)]; end
    tms = [mvmt(ton,2) mvmt(toff,2)]; % this is the on and off times of movement traces under 1.3
    
    thta(:,1) = tX;
    thta(:,2) = data(:,2);
    thta(:,3) = 0;
    
    
    % make a matrix: first column is time, second column is theta integral,
    % third is binary decision of movement below 1.3 threshold.
    for ee = 1:size(tms,1);
        
        p = [];
        p = thta(:,1)>=tms(ee,1) & thta(:,1)<=tms(ee,2) & pts(:,2)>2;
        thta(p,3) = 1;
        
    end
    
    f = thta(:,3) == 1;
    pts(f,2) = 0;
    
    if verbose == 1;
        % find theta above 1.5 and binary 1 for the movement threshold
        targets = [];
        targets = thta(:,2)>1.5 & thta(:,3)==1;
        thta(:,4) = 0;
        thta(targets,4) = 1;
        figure (4)
        plot(tX,thta(:,4)), hold on
        line([tX(1) tX(end)],[1.5 1.5]);
        plot(tX,data(:,2));
        plot(mvmt(:,2),mvmt(:,1),'m')
        set(gca,'ylim',[0 3]);
        title('Movement versus theta power');
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - -
    
    %%%% Now find if theta is low, delta is low, and video is low. This will be
    %%%% quiet wake?
    
    qcombo = []; qt1 = []; qtd = []; qton = []; qtoff = []; qtms = []; qdat = [];
    
    qcombo(:,1) = tX;
    qt1     = mvmt(:,1)>0.3 ;%& mvmt(:,2)<1; % MAKE QT 2 forhigher movement and higher delta
    qtd     = diff(qt1);
    qton    = find(qtd == 1);
    qtoff   = find(qtd == -1);
    if qt1(1) == 1, qton = [1; qton]; end
    if qt1(end) == 1; qtoff = [qtoff; numel(qt1)]; end
    qtms = [mvmt(qton,2) mvmt(qtoff,2)]; % this is the on and off times of movement traces between 0.15 and 0.5
    
    qdat(:,1)  = tX;
    qdat(:,2)  = data(:,2);
    qdat(:,3)  = data(:,1);
    qdat(:,4)  = 0;
    
    
    % make a matrix: first column is time, second column is theta integral,
    % third is binary decision of movement below 1.3 threshold.
    for ee = 1:size(qtms,1);
        
        pq = [];
        pq = qdat(:,1)>=qtms(ee,1) & qdat(:,1)<=qtms(ee,2) & pts(:,2)<3 & qdat(:,2)<1 & qdat(:,3)<1.5;
        qdat(pq,4) = 1;
        
    end
    
    qf = qdat(:,4) == 1;
    pts(qf,2) = 5;
    
    if verbose == 1;
        figure(5)
        plot(tX,qf)
        axis([tX(1) tX(end) 0 4]);
        title('quiet wake identified by all traces being low power...');
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - -
    %%%% SAME THIS AS TWO ITERATIONS ABOVE except now look for medium delta and
    %%%% medium/high movement
    
    mcombo = []; mt1 = []; mtd = []; mton = []; mtoff = []; mtms = []; mdat = [];
    
    mcombo(:,1) = tX;
    mt1     = mvmt(:,1)>1.25 ;%& mvmt(:,2)<1; % MAKE QT 2 forhigher movement and higher delta
    mtd     = diff(mt1);
    mton    = find(mtd == 1);
    mtoff   = find(mtd == -1);
    
    if mt1(1) == 1, mton = [1; mton]; end
    if mt1(end) == 1; mtoff = [mtoff; numel(mt1)]; end
    
    mtms = [mvmt(mton,2) mvmt(mtoff,2)]; % this is the on and off times of movement traces greater than 1.5
    
    mdat(:,1)  = tX;
    mdat(:,2)  = data(:,2);
    mdat(:,3)  = data(:,1);
    mdat(:,4)  = 0;
    
    
    % make a matrix: first column is time, second column is theta integral,
    % third is binary decision of movement below 1.3 threshold.
    for ee = 1:size(mtms,1);
        
        pm = [];
        pm = mdat(:,1)>=mtms(ee,1) & mdat(:,1)<=mtms(ee,2) & pts(:,2)<3 & mdat(:,2)<1.25 & mdat(:,3)<2.5;
        mdat(pm,4) = 1;
        
    end
    
    mf = mdat(:,4) == 1;
    pts(mf,2) = 5;
    
    if verbose == 1;
        figure(5)
        plot(tX,mf)
        axis([tX(1) tX(end) 0 4]);
        title('quiet wake identified by all traces being low power...');
    end
    
    
    wakeID = pts;
else
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % THIS IS CASE IN WHICH THERE IS NO MOVEMENT DATA (i.e. CAMERA IS
    % OFFLINE). HERE, USE THE EMG DATA INSTEAD.
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
    figure (2)
    plot(data(:,3)), hold on
    plot(data(:,1),'r');
    plot(data(:,2),'b');
    hold off
    pts(:,1)    = 1:size(tX,2);
    pts(:,2)    = zeros;
    
    
    EMGAWTHRESH = 1.0;
    

    
    
    % Code active waking based on the EMGAWTHRESH but NOT IF DELTA IS HIGH!
    pts ( data(:,3)>EMGAWTHRESH & data(:,1)<1 & data(:,2)<0.9 ,2 ) = 4; % set to AW
    
    pts( pts(:,2) == 4 & data(:,1)>1.75 & data(:,2)>0.75, 2 ) = 2;
    
    
    
    
    if any( unique(pts(:,2)) == 4 | unique(pts(:,2)) == 5 );
        
        % find the onset and offset of "3" states
        idx3    = pts(:,2) == 3;
        d3      = diff(idx3);
        on3     = find(d3 == 1) + 1;
        off3    = find(d3 == -1) + 1;
        
        % Correct for case in which the block ends with a "3" state
        if size(on3,1) > size(off3,1) && idx3(end) == 1;
            off3 = [off3; numel(idx3)];
        end
        
        
        for ee = 1:size(on3,1);
            
            if on3(ee) - 20 < 0;
                prior = 0;
            else
                prior   = mode( pts(on3(ee)-20 : on3(ee),2) );
            end
            
            if off3(ee) + 20 > numel(pts(:,2));
                post = 0;
            else
                post    = mode( pts(off3(ee) : off3(ee)+20,2) );
            end
            
            pvals = [prior post];
            
            % case in which the "3" is surrounded by nothing.
            if post == 0 && prior == 0;
                pts(on3(ee):off3(ee),2) = 0;
            end
            
            % case in which "3" is flanked by a zero on one side and a "4" or
            % "5" on the other side.
            if any(pvals == 0) && post~=prior
                pts(on3(ee):off3(ee),2) = pvals(pvals~=0);
            end
            
            % case in which "3" is flanked on each side by incongruent non zero elements.
            if all(pvals~=0) && prior~=post;
                halfway = round(on3(ee) + (off3(ee)-1 - on3(ee))/2 );
                pts(on3(ee):halfway)    = prior;
                
                pts(halfway+1:off3(ee),2) = post;
                
            end
            
        end
        
    end
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - -
    
    %%%% Now find if theta is low, delta is low, and video is low. This will be
    %%%% quiet wake?
    pts( data(:,1)>2 & data(:,2)<2 & data(:,3)<1, 2 ) = 5;
    
    
    wakeID = pts;
    
    
    
    
end