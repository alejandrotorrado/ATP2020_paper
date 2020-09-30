function [trips] = find_sleep_triplets(st)

codes = st(:,1);
times = st(:,2);

trips = [];


for ss = 1:numel(codes)-2
    
    s0 = codes(ss);
    
    if s0 == 2 % if find a NREM episode
        
        % check if it is followed by REM and NREM
        next = codes(ss+1);
        nextnext = codes(ss+2);
        
        if next == 1 && nextnext == 2
            
            trips(end+1,1) = ss;
            trips(end,2) = ss+1;
            trips(end,3) = ss+2;
            
        end
            
        
    end
end

if trips(end,3) == size(st,1)
    trips(end,:) = [];
end