function[datback] = flicker(dataidx,minT,time)
%
% [DATBACK] = FLICKER(S_CODE, DATAIDX, TPTS, PTS, MINT)
%
% Iteratively scan though activity/arousal state codes for instances of
% short epochs that represent noise (effective). This will loop through the
% data and try to evaluate which state SHOULD be in place during a (short)
% noise epoch. Epochs shorter than MINT are recognized as erronious. This
% code will only examine one state at a time, but can correct short epochs
% as any of the remaining three states. S_CODE dictates which state is
% subject to editing. DATAIDX is a list of  n sample points by 2 columns.
% The first column is the sample index and the second is the state code.
% TPTS are the unixtime sample times for the data in DATAIDX.
%
%   NOTE: REM = 1, NREM = 2, ACTIVE = 4, QUIET = 5
%
% Output:
%   DATBACK is a n samples x 2 matrix. The first column contains the time
%   indices of sampling. The second column contains the updated
%   "de-flickered" state-codes contained in the second column of DATAIDX.
%
%   10/2015 KBH


count = 0;
ee = 0;

while ee == 0;
    
    trans = []; testD = []; tidx = [];
    testD   = diff([dataidx(1,2) ; dataidx(:,2)]);
    trans   = find(testD ~= 0);
    trans   = [1; trans; numel(time)];
    
    tidx(:,1) = trans(1:end-1);
    tidx(:,2) = trans(2:end);
    tidx(:,3) = time(tidx(:,2)) - time(tidx(:,1));

    if tidx(end,1) == size(dataidx,1);
        tidx(end,:) = [];
    end
    
    epoch = find(tidx(:,3)<=minT,1);
    
    if ~isempty(epoch);
        
        disp(tidx(:,1:2));
        
        if epoch~=size(tidx,1);
            count = count+1;
            disp(['Iteration ' num2str(count)]);
            dataidx(tidx(epoch,1):tidx(epoch,2)-1,2) = dataidx(tidx(epoch,2),2);
        else
            count = count+1;
            disp(['Iteration ' num2str(count)]);
            dataidx(tidx(epoch,1):tidx(epoch,2)-1,2) = dataidx(tidx(epoch,1)-1,2);
        end
        

    else
        disp(tidx(:,1:2))
        ee = 1;
    end

end


datback = dataidx;

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% z = 0; count = 0;
% while z == 0;
%     
%     count = count+1;
%     disp(['Flicker fix pass ' num2str(count)]);
%     
%     STATEidx = dataidx(:,2);
%     STATEidx(STATEidx~=s_code) = 0;% set to binary STATE/not STATE
%     
%     ons     = find(diff(STATEidx) == 1);
%     ons     = ons+1;
%     offs    = find(diff(STATEidx) == -1);
%     
%     
%     if size(ons,1)>size(offs,1);
%         offs = [offs; dataidx(end,1)];
%     elseif size(ons,1)<size(offs,1)
%         ons = [dataidx(1,1); ons];
%     end
%     
%     
%     if any( (tpts(offs)-tpts(ons))<minT );
%         
%         thisSTATE = find( (tpts(offs)-tpts(ons))<minT,1);
%         
%         %check to see if there's STATE immediately before or after this bout
%         
%         if ons(thisSTATE)>60 && offs(thisSTATE)<numel(tpts)-60 && thisSTATE ~= 1;
%             % look back and forward
%             b = tpts(ons(thisSTATE+1)) - tpts(offs(thisSTATE)); % time off before next bout
%             bb = tpts(ons(thisSTATE)) - tpts(offs(thisSTATE-1)); % time since last bout
%             
%             if any([b, bb]<minT);
%                 % change the next or prior bout to CURRENT STATE and then reiterate.
%            
%                 if b < minT;
%                    dataidx( ( offs(thisSTATE)+1:ons(thisSTATE+1)-1), 2) = s_code; 
%                 else
%                     keyboard
%                     %check this line of code:
%                     dataidx( ( offs(thisSTATE-1)+1:ons(thisSTATE)-1), 2) = s_code; 
%                 end
%                 
%             else
%                 % figure out which flanking state is most similar to the
%                 % error state and merge the two:
%                 cprior  = dataidx(ons(thisSTATE)-1,2);
%                 cpost   = dataidx(offs(thisSTATE)+1,2);
%                 temp(:,1) = [cprior; cpost];
%                 c = repmat(dataidx(ons(thisSTATE),2),2,1);
%                 temp(:,2) = c;
%                 temp(:,3) = abs(temp(:,1) - temp(:,2) );
%                 dataidx(ons(thisSTATE):offs(thisSTATE),2) = unique(temp(  temp(:,3)==min(temp(:,3)), 1));
%                 
%             end
%             
%             
%             
%         elseif ons(thisSTATE)<60 | thisSTATE == 1;
%             % only look forwards
%             b = tpts(ons(thisSTATE+1)) - tpts(offs(thisSTATE)); % time off before next bout
%             
%             if b<minT;
%                 % change the next bout to CURRENT STATE and then reiterate.
%                 dataidx(offs(thisSTATE)+1:ons(thisSTATE+1)-1,2)   = s_code;
%             else
%                 dataidx(ons(thisSTATE):offs(thisSTATE),2)   = dataidx(offs(thisSTATE)+1,2);
%             end
%             
%         elseif ons(thisSTATE)>60 && offs(thisSTATE)>numel(tpts)-60;
%             % only look backwards
%             keyboard
%             bb = tpts(ons(thisSTATE)) - tpts(offs(thisSTATE-1)); % time since last bout
%             
%         end
%         
%     else
%         % everything has been worked out
%         z = 1;
%         
%     end
%     
% end
% 
% datback = dataidx;