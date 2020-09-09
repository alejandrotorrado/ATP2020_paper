
% FIRST: load prior training dataset, then add a column of zeros to the
% end, and move the quality to that column, then reset the second to last
% as all zeros in order to set up a new column for another factor.



% now go through and load the channel data from each animal contributing to
% p and run this code to then find the right row (that matches the given
% cluster) and add a new factor to it. 
sDir = uigetdir;
%chanfiles = dir([sDir filesep '*channel_*']);
cellfile = dir([sDir filesep '*_MASTER_CELL.mat']);
load([sDir filesep cellfile.name]);
pfile = dir([sDir filesep '*_sorting_p.mat']);
load([sDir filesep pfile.name]);

% load the model
load('/Users/khengen/Google_Drive/Matlab_scripts_11_01_2012/Cell_Quality/RANDOMFOREST_TRAINED/rForest_Mdl_KH72_73_75_SHANK01_AT14_trained.mat')

for ee = 1:size(CELL,2)
    
    
    X = diff(CELL(ee).time);
    h = histc(X,[0:0.001:0.5]);
    
    if sum(h)>=1000
        tmp = gamfit(1:100,[],[],h(1:100)');
        
        target = find(p_save(:,1) == tmp(1) & p_save(:,2) == tmp(2));
        if numel(target)~=1
            disp('ISSUES');
            keyboard
        end
        
        Y = p_save(target,1:14);
        [label,scoretmp,~] = predict(Mdl,Y); % use model to score
        label = str2double (label{1});
        
        OGq = CELL(ee).quality;
        
        CELL(ee).quality = label;
        
        disp(['OG ' num2str(OGq) ' new ' num2str(CELL(ee).quality)]);
        
        score = [];
        score(1,1) = label;
        score(1,2) = scoretmp(label);
        scoretmp2 = scoretmp;
        scoretmp2(label) = -100;
        [scoremax,~] = max(scoretmp2);
        [~,scoremaxlabel] = find(scoretmp2 == scoremax);
        
        if numel(scoremaxlabel)==1
            score(1,3) = scoremaxlabel;
            score(1,4) = scoremax;
        else
            score(1,3) = 0;
            score(1,4) = 0;
        end
        
        
        CELL(ee).score = [];
        CELL(ee).score = score(1,:);
        
    else
        %do nothing
    end
    
end


% 
% count = 0;
% for ii =  1:size(chanfiles,1)
%     
%     if exist([sDir filesep chanfiles(ii).name],'file') == 2
%         
%         load([sDir filesep chanfiles(ii).name]);
%         
%         for ee = 1:size(chdata.block{1}.clust,2)
%        
%             dat = [];
%             dat = chdata.block{1}.clust(ee).time;
%             X = diff(dat);
%             h = histc(X,[0:0.001:0.5]);
%             
%             if sum(h)<1000
%                 % do nothing
%             else
%                 
%                 tmp = gamfit(1:100,[],[],h(1:100)');
%                 
%                 target = find(p(:,1) == tmp(1) & p(:,2) == tmp(2));
%                 if numel(target)~=1
%                     disp('ISSUES');
%                     keyboard
%                 end
%                 
%                 % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%                 % DO THIS FOR UPDATING:
%                 % - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%                  Fs = 1000;           % Sampling frequency
%                 L = 500;             % Length of signal
%                 fhz = Fs*(0:(L/2))/L;
%                 yfft = fft(h(1:end-1));
%                 P2 = abs(yfft/L);
%                 P1 = P2(1:L/2+1);
%                 P1(2:end-1) = 2*P1(2:end-1);
%                 
%                 P1norm = P1/sum(P1);
%                 %plot(fhz,P1norm)
%                 peaktemplate = [0.3868 0.4064 0.4220 0.4371 0.4596 0.4859 0.5210... 
%                     0.5651 0.6203 0.6875 0.7944 1.0000 0.7902 0.6830 0.6138 ...
%                     0.5590 0.5119 0.4745 0.4461 0.4201 0.4005]';
% 
%                 ac_noise  = corr(peaktemplate, P1(20:40)/max(P1(20:40)) );
%                 % - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%   
%                 p(target,14) = ac_noise;
%                 disp(target);
%             end
%             
%         end
% 
%     end
%     
%     
% end
% 
% % Now do something about them!
% fixers = find(p(:,14)>0.5 & p(:,15)<3)
