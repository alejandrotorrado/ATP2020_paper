function [bac] = overclusteringcheck(clusters)
%    OVERCLUSTERINGCHECK - Routine to compare the mean square error of
%    every possible combination of unique cluster average waveforms. If
%    the MSE of a comparison (or set of comparisons) is very small, the
%    two clusters will be merged into a single cluster. This algorithm
%    has the capacity to detect elements shared between more than one
%    pair and merge the unique elements of the multiple pairs into a
%    group.
%
%    INPUTS
%        clusters - the data contained in "channel(cc).block{bb}.clust"
%        which should have the following fields: idx, points, ctr,
%        dists, clustvar, samples, meantrace, tracestdev, and time. This
%        structure is output from the file klustakwikRJ.m, and perhaps
%        others. Reference those scripts for a more complete description
%        of the contents of the fields of "clusters".
%
%    OUTPUTS
%        bac - either a modified or unchanged version of clusters.
%        Changes reflect the merging of highly similar clusters input in
%        "clusters" and deletion of the contributing clusters following
%        the merge.
%

nsampl  = size(clusters(1).meantrace,2);
meanwfs = cell2mat({clusters.meantrace});
meanwfs = reshape(meanwfs,nsampl,[]);

% get all possible combinations of the clusters:
c       = combnk(1:size(clusters,2),2);
c(:,3)  = zeros;

% calculate the mean squared error for each possible pair of clusters
for ee = 1:size(c,1);
    c(ee,3) = meansquareerr( meanwfs(:,c(ee,1)), meanwfs(:,c(ee,2)) );
end

% find rows (pairs) that have MSE<20;
d       =  c(c(:,3)<20,1:2);
verbose = 0;

if ~isempty(d);
    disp(['Overclustering detected. Clusters involved: ' num2str(unique(d(:))') '.'] );
    
    if verbose;
        plotty = unique(d(:));
        for pp = 1:size(plotty,1);
            figure(888); hold on
            plot(clusters(plotty(pp)).meantrace);
        end
        title('Mean traces of merging clusts');
        hold off
    end
    
    % make each row into a cell
    temp    = mat2cell(d,ones(size(d,1),1));
    
    % look for intersection between rows. if it's there, merge them.
    rows    = combnk(1:size(temp,1),2);
    
    temp2 = {}; count = 0; temptracker = zeros(size(temp,1),1);
    for ee = 1:size(rows,1);
        
        together = intersect(temp{rows(ee,1)}, temp{rows(ee,2)});
        
        if ~isempty(together);
            temptracker(rows(ee,1)) = 1;
            temptracker(rows(ee,2)) = 1;
            if isempty(temp2)
                count = count+1;
                temp2{count} = unique([temp{rows(ee,1)} temp{rows(ee,2)}]);
            else
                % look to see if the intersected pair fits with a previous set
                % of matched pairs. if so, merge them. if not, create a new
                % entry in temp2
                ii = 0; icount = 0; hook = 0;
                while ii == 0;
                    icount = icount+1;
                    catcher = intersect( temp2{icount}, unique([temp{rows(ee,1)} temp{rows(ee,2)}]) );
                    if ~isempty(catcher);
                        temp2{icount} = [temp2{icount} unique([temp{rows(ee,1)} temp{rows(ee,2)}])];
                        ii = 1;
                        hook = 1;
                    end
                    
                    if icount == size(temp2,2);
                        ii = 1;
                    end
                end
                
                if hook == 0;
                    count = count+1;
                    temp2{count} = unique([temp{rows(ee,1)} temp{rows(ee,2)}]);
                end
            end
        end
    end
    
    % find any pairs that had no shared values with other rows
    if any(temptracker == 0);
        
        adds = find(temptracker == 0);
        
        for ee = adds';
            count = count+1;
            temp2{count} = temp{ee};
        end
    end
    
    
    for ee = 1:size(temp2,2);
        
        d2{ee} = unique(temp2{ee});
        
    end
    
    % Build the new cluster structures and then replace the contributing
    % clusters in the dataset.
    
    for cc = 1:length(d2)
        clear newclust
        clulist = d2{cc};
        
        % build a matrix of whose mean trace represents the mean WF of the new
        % cluster
        matmatrix = [];
        for ee = 1:length(clulist);
            matrix      = [];
            matrix      = repmat(clusters(clulist(ee)).meantrace,size( clusters(clulist(ee)).idx,1 ),1 );
            matmatrix   = [matmatrix ; matrix];
        end
        
        totaln  = size(matmatrix,1);
        
        combvar = sum(var(cell2mat({clusters(clulist).points}')));

        newclust.idx        = cell2mat({clusters(clulist).idx}');
        newclust.points     = cell2mat({clusters(clulist).points}');
        newclust.ctr        = mean( newclust.points );
        newclust.dists      = mahaldist(newclust.points);
        newclust.clustvar   = combvar;
        newclust.samples    = 'merged'; % since you don't have the OG samples, this is a pointless field.
        newclust.meantrace  = mean(matmatrix);
        newclust.tracestdev = NaN(size(newclust.meantrace));
        newclust.time       = cell2mat({clusters(clulist).time}');
        if isfield(clusters(clulist),'channel')
            newclust.channel    = unique(cat(1,clusters(clulist).channel));
        end
        newclust.quality    = [];
        try
            clusters(end +1)    = newclust;
        catch
            keyboard
        end
    end
   
    clusters(cell2mat(d2)) = [];
    
else
    disp('No overclustering detected. Leaving current clusters intact.');
end

bac = clusters;



