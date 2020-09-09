function clust = klustakwikRJ_hpc( data, minclusters, maxclusters, displayflag, fileNo )
global myname;

% data argument is the PCA output

% optional arguments
if nargin == 4
    fileNo = '0';
end

% set up file names and parameters
dim = size(data, 2);
basename    = ['KlustTempData_' fileNo];
ElecNo      = 1;                            % run identifier, any integer, but required by KlustaKwik
fetname     = [basename '.fet.' num2str(ElecNo) ];
cluname     = [basename '.clu.' num2str(ElecNo) ];


% IMPORTANT
% the SubSetAt parameter is key. Subsetting means using a subset of the
% spikes on a channel to cluster the channel. KlustaKwik then assigns the
% rest of the data to the clusters based on similarity. This speeds up the
% process - it would be impossible to cluster 200+ hours of recordings and
% 300 million spikes at once. Thus, only up to 5 million spikes per channel
% are used.
SubSetAt = 5e6; % number of spikes to use for clustering
nspikesTotal = size(data,1);
% but do this only if there are too many spikes
if nspikesTotal > SubSetAt
    SubSetRate = round(nspikesTotal/SubSetAt);
    fprintf('Subsetting at rate of 1/%d\n', SubSetRate);
else
    SubSetRate = 1;
end

% write pca data to features file
formatstr = [ repmat( ['%2.2f\t'], 1, dim ) '\n' ]; 
fid = fopen(fetname,'wt','b');
fprintf(fid,'%u\n', dim);
fprintf(fid, formatstr, data');
fclose(fid);

% Make a shell call to KlustaKwik
if strcmp(displayflag, 'verbose') || strcmp(displayflag, 'both') 
    commandstr = [ ... % sdv modified
                    char(32) basename char(32) num2str(ElecNo) ... 
                    ' -UseDistributional '      '0'...
                    ' -MinClusters '            num2str(minclusters) ... 
                    ' -MaxClusters '            num2str(maxclusters) ... 
                    ' -MaxPossibleClusters '    num2str(maxclusters) ... 
                    ' -PenaltyK '               '0'...
                    ' -PenaltyKLogN '           '1'...
                    ' -Subset '                  num2str(SubSetRate)...
                    ' -Screen '                 '1';
                   ] ;      % outputting text to command window - the log is saved as *.klg
else
    commandstr = [ ... % sdv modified
                    char(32) basename char(32) num2str(ElecNo) ... 
                    ' -UseDistributional '      '0'...
                    ' -MinClusters '            num2str(minclusters) ... 
                    ' -MaxClusters '            num2str(maxclusters) ... 
                    ' -MaxPossibleClusters '    num2str(maxclusters) ... 
                    ' -PenaltyK '               '0'...
                    ' -PenaltyKLogN '           '1'...
                    ' -Subset '                  num2str(SubSetRate)...
                    ' -Screen '                 '0'...
                    %' > /dev/null';
                    ] ;      % avoid outputting text to command window - the log is saved as *.klg
end

% run klustakwik using wrapper
klustakwik_matlab_hpc(commandstr);

% open the output cluster file
fid = fopen(cluname, 'r');

% read output
if fid<0, error(['Error opening file ' cluname '.']); end;
   clustnum    = fread(fid, inf, '*char');

clustnum   = strsplit(clustnum');
clustnum   = str2double(clustnum');

% First entry is number of clusters found - cluster number 1 is noise or outliers, and may be empty
nclusts     = clustnum(1);
clustnum(1) = [];                               % Now delete that first entry so all remaining entries are real cluster numbers
fclose(fid);

disp([num2str(nclusts) ' clusters found.']); 

% format the KlustaKwik output for our purposes
for c = 1:nclusts
    clust(c).idx        = find( clustnum == c );
    clust(c).idx        = clust(c).idx(:);

    if ~isempty( clust(c).idx )
        %         Use raw data points
        clust(c).points     = data( clust(c).idx, 1:3 );
        clust(c).ctr        = mean( clust(c).points );
        clust(c).dists      = mahaldist(clust(c).points);
        clust(c).clustvar   = sum(var(clust(c).points));
        
        % legacy display flag - DO NOT USE ON CLUSTER
        if strcmp(displayflag, 'display') | strcmp(displayflag, 'both')
            %             clust(c).plothandle = plot3( clust(c).points(:, 1) , clust(c).points(:, 2) , clust(c).points(:, 3), '.' , 'color', colordc(c,:)); hold on
            clust(c).plothandle = plot3( P.Proj(clust(c).idx, 1) , P.Proj(clust(c).idx, 2) , P.Proj(clust(c).idx, 3), '.' , 'color', colordc(c,:)); hold on
            axis vis3d; rotate3d on; grid on
            drawnow
        end
    end
end







 
