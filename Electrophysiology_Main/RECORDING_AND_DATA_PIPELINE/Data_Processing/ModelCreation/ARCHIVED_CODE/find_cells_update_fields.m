
% FIRST: load prior training dataset, then add a column of zeros to the
% end, and move the quality to that column, then reset the second to last
% as all zeros in order to set up a new column for another factor.



% now go through and load the channel data from each animal contributing to
% p and run this code to then find the right row (that matches the given
% cluster) and add a new factor to it. 
sDir = uigetdir;
chanfiles = dir([sDir filesep '*channel_*']);

count = 0;
for ii =  1:size(chanfiles,1)
    
    if exist([sDir filesep chanfiles(ii).name],'file') == 2
        
        load([sDir filesep chanfiles(ii).name]);
        
        for ee = 1:size(chdata.block{1}.clust,2)
       
            dat = [];
            dat = chdata.block{1}.clust(ee).time;
            X = diff(dat);
            h = histc(X,[0:0.001:0.5]);
            
            if sum(h)<1000
                % do nothing
            else
                
                temp = gamfit(1:100,[],[],h(1:100)');
                
                target = find(p(:,1) == temp(1) & p(:,2) == temp(2));
                if numel(target)~=1
                    disp('ISSUES');
                    keyboard
                end
                
                % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                % DO THIS FOR UPDATING:
                % - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                 Fs = 1000;           % Sampling frequency
                L = 500;             % Length of signal
                fhz = Fs*(0:(L/2))/L;
                yfft = fft(h(1:end-1));
                P2 = abs(yfft/L);
                P1 = P2(1:L/2+1);
                P1(2:end-1) = 2*P1(2:end-1);
                
                P1norm = P1/sum(P1);
                %plot(fhz,P1norm)
                peaktemplate = [0.3868 0.4064 0.4220 0.4371 0.4596 0.4859 0.5210... 
                    0.5651 0.6203 0.6875 0.7944 1.0000 0.7902 0.6830 0.6138 ...
                    0.5590 0.5119 0.4745 0.4461 0.4201 0.4005]';

                ac_noise  = corr(peaktemplate, P1(20:40)/max(P1(20:40)) );
                % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
                p(target,14) = ac_noise;
                disp(target);
            end
            
        end

    end
    
    
end

% Now do something about them!
fixers = find(p(:,14)>0.5 & p(:,15)<3)
