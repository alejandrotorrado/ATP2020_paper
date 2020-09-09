function  [datback] =  REMQWfilter(~, hseq, d3, pX, verbose)




statetimes = d3;

% find state transitions:
testD   = diff(statetimes(:,2));
killem  = find(testD == 0);
killem  = killem+1;
statetimes (killem,:) = [];

r = find(statetimes(:,2) == 1);Ron = statetimes(r,1);


if ~isempty(r); % make sure there are REM bouts in this block of data
    if Ron(end) == statetimes(end,1);
        
        Roff= statetimes(r(1:end-1)+1 ,1);
        Roff = [Roff; d3(end,1)];
    else
        Roff= statetimes(r+1 ,1);
    end
    
    if verbose == 1;
        figure(55555);
        plot(pX,hseq(:,1),'color',[0 0.4 1]), hold on
        plot(pX,hseq(:,2),'b');
        plot(pX,hseq(:,3),'color',[1 0.4 0])
        legend('Delta','Theta','EMG');
        
        for ee = 1:size(Roff,1);
            h(ee) = line([pX(Ron(ee)) pX(Ron(ee))],[0 4]);
            set(h(ee),'color','g');
            g(ee) = line([pX(Roff(ee)) pX(Roff(ee))],[0 4]);
            set(g(ee),'color','r');
        end
        hold off
    end
    
    %calculate REM theta powers
    for ee = 1:size(Roff,1);
        tpow(ee) = mean(hseq(Ron(ee):Roff(ee),2));
        
        bins = [];
        bins = Ron(ee):100: Roff(ee);% Bin at apprx 50 sec
        
        tpow2 = []; dpow = [];
        for ii = 1:size(bins,2)-1;
            
            tpow2(ii)    = mean(hseq(bins(ii):bins(ii+1),2) );
            dpow(ii)     = mean(hseq(bins(ii):bins(ii+1),1) );
            
            
        end
        
        
        % look for chunks of REM that have really low delta and theta and
        % reclassify these as QW
        
        test = tpow2<1 & dpow<1.5;
        
        t1 = diff(test);
        
        % - ONLY USE THOSE THAT HAVE MORE THAN 2 CONTIGUOUS BLOCKS
        oo = [];
        if sum(tpow2<1 & dpow<1.5) >2;
            
            oo(:,1) = bins(tpow2<1 & dpow<1.5);
            oo(:,2) = bins(tpow2<1 & dpow<1.5) + 100;
            
            for yy = 1:size(oo,1);
                d3(oo(yy,1):oo(yy,2),2) = 5; % convert this to quiet wake
            end
        end
        
    end
    
end

datback = d3;





