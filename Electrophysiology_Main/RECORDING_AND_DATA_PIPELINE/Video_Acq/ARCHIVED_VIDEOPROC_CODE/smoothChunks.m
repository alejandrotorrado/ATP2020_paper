function [smoothed_mvt] = smoothChunks(mvt,sm_factor,chunk_size)

rows = size(mvt,1);

n_chunks = ceil(rows/chunk_size);

for chunk = 1:n_chunks
    fprintf('Smoothing chunk %u out of %u.\n', chunk, n_chunks);
    
    start_pos   = chunk_size*(chunk-1)+1;
    end_pos     = chunk_size*chunk;
    if end_pos > rows, end_pos = rows; end
    
    mvt_chunk = mvt(start_pos:end_pos,1);
    
    t0 = tic;
    smooth_chunk = smooth(mvt_chunk,sm_factor,'rloess');
    t1 = toc(t0);
    fprintf('That took %.2f seconds.\n\n',t1);

    smooth_chunk(smooth_chunk < 0) = 0;
    
    smoothed_mvt(start_pos:end_pos,1) = smooth_chunk;
    
    clear mvt_chunk smooth_chunk start_pos end_pos
end




