function [datback] = alignandcluster (sweepdat,thischan,yy,houredges,subidx,scrubbedidx,timestamps,doInterp,interp_samp)

% interpolate and re-center the spikes:
if doInterp
    fprintf('\nDoing interpolation. Current number of samples in spikes: %u.\n',size(sweepdat,2));
    
    tic
    x   = 1:size(sweepdat,2);
    xq  = 0: 1/3 :size(sweepdat,2);
    disp('Set up finished. Commencing interpolation.')
    sweepdat = double(sweepdat);
    try
        sweepdat = interp1(x,sweepdat',xq,'spline');
    catch
        error('alignandcluster failed at sweepdat inteprolation');
    end
    sweepdat = sweepdat';
    
    
    
    [~,minidx2]     = min(sweepdat(:,10:30),[],2); % this is 20:50 if you're interpolating at 5x. 
    minidx2         = minidx2 + 9;
    swprng          = -18:40;
    
    
    sweepdat        = alignrows(sweepdat, minidx2, swprng, size(sweepdat,2)-1);
    
    toc
else
    % align
    if interp_samp == 97
        [~,minidx2]     = min(sweepdat(:,10:30),[],2); % this is 20:50 if you're interpolating at 5x.
        minidx2         = minidx2 + 9;
        swprng          = -18:40;
    elseif interp_samp == 161
        [~,minidx2]     = min(sweepdat(:,20:50),[],2); % this is 20:50 if you're interpolating at 5x.
        minidx2         = minidx2 + 19;
        swprng          = -30:60;
    end
    
    
    sweepdat        = alignrows(sweepdat, minidx2, swprng, size(sweepdat,2)-1);
end
%  - - - - - - - - - - - - - do PCA - - - - - - - - - - - -
% newindx         = origindx;
%     Do the PCA
try
    % fprintf('Doing PCA on sweepdat. Size: %u,%u.\n',size(sweepdat,1),size(sweepdat,2));
    fprintf('Doing PCA.\n');
    pcat0 = tic;
    P = pca_hpc( sweepdat );         % Principal Component Analysis
    % disp(size(P.V));
    pcat1 = toc(pcat0);
    fprintf('PCA took %.2f seconds.\n\n',pcat1);
catch
    disp('PCA issue');
end

% fprintf('Got P. Size: %u,%u.\n',size(P,1),size(P,2));

for n = 1:size(P.V)     % Adjust Eigenvectors so maximum deflection is always downward (to aid visual comparison)
    mx = max(P.V(:,n));
    mn = min(P.V(:,n));
    if mx > abs(mn); P.V(:,n) = -P.V(:,n); end
end

minclusts   = 1;
if isempty(houredges)
    maxclusts = 9;
else
    maxclusts = 8;
end
ClusterSet  = P.Proj(:, 1:4);%KBH changed this from 3 to 4 on 4/26/12

chanblockcode = [num2str(thischan) '_' num2str(yy)];

fprintf('Starting clustering!\n');
block.clust     = klustakwikRJ_hpc( ClusterSet, minclusts, maxclusts, 'notverbose', chanblockcode );


for gg = 1:size(block.clust,2);
    
    if size(block.clust(gg).idx,1)>100;
        % pull 100 random waveforms and save them
        y = datasample(block.clust(gg).idx,100,'Replace',false);
        block.clust(gg).samples = sweepdat(y,:);
    else
        y = block.clust(gg).idx;
        block.clust(gg).samples = sweepdat(y,:);
    end
    
    block.clust(gg).meantrace   = mean(sweepdat(block.clust(gg).idx,:));
    block.clust(gg).tracestdev  = std( sweepdat( block.clust(gg).idx,: ) );
    
    % update idx to correspond to the original spikes files:
    if ~isempty(subidx)
        block.clust(gg).idx     = scrubbedidx( subidx( block.clust(gg).idx ) ); % this is the index of spikes corresponding to the binary files at the start
    else
        block.clust(gg).idx     = scrubbedidx( block.clust(gg).idx ) ; % indices of spikes corresponding to binary file
    end
    
    block.clust(gg).time    = timestamps( block.clust(gg).idx );
    
end

% Add more fields to the block variable for later manipulation:
if isempty(houredges)
    block.start = 0;
    block.end = timestamps(end);
else
    block.start = houredges(yy);
    block.end   = houredges(yy+1);
end

block.interpSamples = interp_samp;
datback = block;
