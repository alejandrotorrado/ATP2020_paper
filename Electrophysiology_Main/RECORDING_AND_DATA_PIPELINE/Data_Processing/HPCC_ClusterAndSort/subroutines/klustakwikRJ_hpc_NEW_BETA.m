function clust = klustakwikRJ_hpc_NEW_BETA( data, minclusters, maxclusters, displayflag, fileNo, subset_n )
global myname;
%cd('C:\Users\khengen.USERS\Documents\MATLAB\');
if nargin == 4
    fileNo = '0';
    subset_n = 1e6;
elseif nargin == 5
    subset_n = 1e6;
end


% cluster data is up sampled at 5x

dim = size(data, 2);
basename    = ['KlustTempData_' fileNo];
ElecNo      = 1;                            % run identifier, any integer, but required by KlustaKwik
fetname     = [basename '.fet.' num2str(ElecNo) ];
cluname     = [basename '.clu.' num2str(ElecNo) ];


% pick .dat file for the channel and check nspikes
SubSetAt = subset_n;
nspikesTotal = size(data,1);
if nspikesTotal > SubSetAt;
    SubSetRate = round(nspikesTotal/SubSetAt);
    fprintf('Subsetting at rate of 1/%d\n', SubSetRate);
else
    SubSetRate = 1;
end
 
% % ADD THIS AS AN INPUT VARIABLE
% try
% if size(data,1)>1e6;
%    %factor   = (size(data,1)/1e6);
%    %subset   = round(rand( round(size(data,1)/factor),1)*size(data,1));
%    subset   = randsample([1:size(data,1)],1e6);
%    OGdata   = data;
%    data     = data(subset,:);
% end
% catch
%     disp('I''m stuck on data subsetting');
%     keyboard
% end


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


klustakwik_matlab_hpc(commandstr);


fid = fopen(cluname, 'r');

if fid<0, error(['Error opening file ' cluname '.']); end;
   clustnum    = fread(fid, inf, '*char');

clustnum   = strsplit(clustnum');
clustnum   = str2double(clustnum');
%    clustnum    = str2num( clustnum(1:1:end) );     % skip the linefeed characters ***This is toms Change***
%         clustnum    = str2num( clustnum(1:2:end) );     % skip the linefeed characters

% First entry is number of clusters found - cluster number 1 is noise or outliers, and may be empty
nclusts     = clustnum(1);
clustnum(1) = [];                               % Now delete that first entry so all remaining entries are real cluster numbers
fclose(fid);

disp([num2str(nclusts) ' clusters found.']);


for c = 1:nclusts
    clust(c).idx        = find( clustnum == c );
    clust(c).idx        = clust(c).idx(:);

    if ~isempty( clust(c).idx )
        %         Use raw data points
        clust(c).points     = data( clust(c).idx, 1:3 );
        clust(c).ctr        = mean( clust(c).points );
        clust(c).dists      = mahaldist(clust(c).points);
        clust(c).clustvar   = sum(var(clust(c).points));
        
        % isolation distance
        clust(c).IsoDist    = Iso_Distance(data, clust(c).idx);
        
        % L-ratio
        LR_temp             = calc_L_Ratio(data, clust(c).idx);
        clust(c).L_ratio    = LR_temp.Lratio;
        
        if strcmp(displayflag, 'display') | strcmp(displayflag, 'both')
            %             clust(c).plothandle = plot3( clust(c).points(:, 1) , clust(c).points(:, 2) , clust(c).points(:, 3), '.' , 'color', colordc(c,:)); hold on
            clust(c).plothandle = plot3( P.Proj(clust(c).idx, 1) , P.Proj(clust(c).idx, 2) , P.Proj(clust(c).idx, 3), '.' , 'color', colordc(c,:)); hold on
            axis vis3d; rotate3d on; grid on
            drawnow
        end
    end
end







 
