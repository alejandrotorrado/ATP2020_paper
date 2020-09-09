function [datback, p_out, trimvec] = alignandcluster_HPC (sweepdat,thischan,yy,houredges,subidx,scrubbedidx,timestamps,doInterp,interp_samp, Mdl, nminrng, nswprng,maxclusts)

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
        fprintf('sweepdat keyboard');
        keyboard;
    end
    sweepdat = sweepdat';
    
    
    
    [~,minidx2]     = min(sweepdat(:,nminrng),[],2); % this is 20:50 if you're interpolating at 5x.
    minidx2         = minidx2 + nminrng(1) - 1;
    
    sweepdat        = alignrows(sweepdat, minidx2, nswprng, size(sweepdat,2)-1);
    
    toc
else
    % align
	[~,minidx2] = min(sweepdat(:,nminrng),[],2);
	minidx2    = minidx2 + nminrng(1) - 1;

    %{
    if interp_samp == 97
        [~,minidx2]     = min(sweepdat(:,10:30),[],2); % this is 20:50 if you're interpolating at 5x.
        minidx2         = minidx2 + 9;
        swprng          = -18:40;
    elseif interp_samp == 161
        [~,minidx2]     = min(sweepdat(:,20:50),[],2); % this is 20:50 if you're interpolating at 5x.
        minidx2         = minidx2 + 19;
        swprng          = -30:60;
    end
    %}
    
    sweepdat        = alignrows(sweepdat, minidx2, nswprng, size(sweepdat,2)-1);
end

%  - - - - - - - - - - - - - do PCA - - - - - - - - - - - -
% newindx         = origindx;
%     Do the PCA
try
    disp('Doing PCA.');
    tic
    P = pca_hpc( sweepdat );         % Principal Component Analysis
    toc
catch
    disp('PCA issue');
    keyboard;
end

for n = 1:size(P.V)     % Adjust Eigenvectors so maximum deflection is always downward (to aid visual comparison)
    mx = max(P.V(:,n));
    mn = min(P.V(:,n));
    if mx > abs(mn); P.V(:,n) = -P.V(:,n); end
end

minclusts   = 1;
% if isempty(houredges)
%     maxclusts = 9;
% else
%     maxclusts = 8;
% end
ClusterSet  = P.Proj(:, 1:4);%KBH changed this from 3 to 4 on 4/26/12

chanblockcode = [num2str(thischan) '_' num2str(yy)];

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
try
    % Run supervised training autosorting here and determine whether a cluster
    % needs trimming:
    fprintf('Performing automated cluster evaluation.');
    plotty = 0; % flag to plot or not
    [newblock, p_out, trimvec] = trained_auto_sorting_ncol20(Mdl,block,interp_samp,interp_samp,nswprng,nminrng,plotty);
    
catch
	disp('Error in alignandcluster.m');
    keyboard
end
% Add more fields to the block variable for later manipulation:
if isempty(houredges)
    newblock.start = 0;
    newblock.end = timestamps(end);
else
    newblock.start = houredges(yy);
    newblock.end   = houredges(yy+1);
end

newblock.interpSamples = interp_samp;
datback = newblock;
