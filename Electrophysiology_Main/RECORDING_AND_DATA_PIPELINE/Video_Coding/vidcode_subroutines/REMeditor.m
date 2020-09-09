function [datout] = REMeditor (statecodes, signals, verbose, times, mvmt)
% REM sleep is easy to confuse for waking states. Try to find inappropriate
% occurences of REM here, and then iron them out


b = statecodes(:,2) == 1;
if b(end) == 1; b(end) = 0; end
c = diff(b);
rem_on = find(c == 1);
rem_off = find(c == -1);
if statecodes(1,2) == 1; rem_on = [1; rem_on]; end


% look for flanking wake states:

checkRon = rem_on(rem_on~=1);
for ee = 1:length(checkRon);
    
    prior = statecodes(checkRon(ee)-1,2);
    
    if prior>3;
        disp('caught example');
        
        thisOff = rem_off(find(rem_off>checkRon(ee),1));
        
        post = statecodes(thisOff+1,2);
        
        if verbose == 1;
            figure(999)
            plot(statecodes(:,1),signals(:,:)), hold on
            g = line([checkRon(ee) checkRon(ee)],[0 4]);
            set(g,'color','g');
            h = line([thisOff thisOff],[0 4]);
            set(h,'color','m');
            axis([statecodes(1,1) statecodes(end,1) 0 4]);
            legend('Delta','Theta','EMG')
            hold off
        end
        
        thesePoints = (checkRon(ee):thisOff)';
        checkEMG = mean(signals(thesePoints,3));
        
        if post > 3;
            
            if checkEMG >1.5;
                statecodes(thesePoints,2) = 4;
            else
                
                statecodes(thesePoints(1: round(0.5*size(thesePoints,2))),2) = prior;
                statecodes(thesePoints(round(0.5*size(thesePoints,2))) : thisOff,2) = post;
                
            end
            
        else
            
            checkdelta = mean(signals(thesePoints,1));
            
            if checkdelta<2 && checkEMG<1.25; % this will be quiet waking
                
                statecodes(thesePoints,2) = 5;
                
            elseif checkdelta>=2 && checkEMG<1; % this will be NREM
                
                statecodes(thesePoints,2) = 2;
                
            elseif checkdelta>2 && checkEMG>1;
                
                % check the movement during this period. it's note REM
                % because the preceeding state was waking, but it's unclear
                % if it the animal was in NREM or quiet waking at this
                % point. Use video tracking to disambiguate this.
                test = mean( mvmt( mvmt(:,2)>= times(checkRon(ee)) & mvmt(:,2)<=times(thisOff),1) );
                if test>0.5;
                    statecodes(thesePoints,2) = 5; % quiet waking
                else
                    statecodes(thesePoints,2) = 2; % NREM
                end
                
                
                
            end
            
            %             conNREM = signals(thesePoints ,1)>=1.5;
            %             statecodes(thesePoints(conNREM),2) = 2;
            
        end
        
    end
    
    
    
    
end

datout = statecodes;
