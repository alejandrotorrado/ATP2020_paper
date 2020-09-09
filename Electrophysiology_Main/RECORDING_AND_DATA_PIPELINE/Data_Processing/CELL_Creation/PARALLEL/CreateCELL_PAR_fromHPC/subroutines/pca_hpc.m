function P = pca_hpc( data )

% P = PCA( data );
% 
% Perform Principal Component Analysis.
% DATA is an NxM matrix of N observations (e.g., traces) and M variables (e.g., sample points).
% P is a struct containing:
%     C, the covariance matrix
%     V, the eigenvector matrix
%     D, the eigenvalue matrix
%     Proj, the projection of all original data points into the eigenvector space
%     
% M Jones, Aug 2008
%     
% example:
% dt = 1000e-6;                   % dt = 100 us (sampling at 10 kHz)
% t = [0:dt:0.1];                 % each sweep = 0 to 100 msec;
% f1 = 13;                        % f1 = 1 Hz % Three Underlying Waveforms
% f2 = 23;                        % f2 = 2 Hz
% f3 = 37;                        % f3 = 3 Hz
% s1 = sin(f1.*t.*2.*pi) .* exp(-f1.*t); 
% s2 = sin(f2.*t.*2.*pi) .* exp(-f2.*t); 
% s3 = sin(f3.*t.*2.*pi) .* exp(-f3.*t); 
% D1 = [rand(33,1)] * s1;                  % Create multiple waveforms add noise
% D2 = [rand(33,1)] * s2;
% D3 = [rand(33,1)] * s3;
% D1 = D1 + 0.1.*randn(size(D1));
% D2 = D2 + 0.1.*randn(size(D2));
% D3 = D3 + 0.1.*randn(size(D3));
% D = [D1; D2; D3];                         % concatenate
% D = D(randperm(size(D,1)), :);          % scramble
% P = pca( D );
% figure
% subplot(2, 2, 1); plot( t, D', '-' );                                                title('Raw Data');              axis square
% subplot(2, 2, 2); imagesc( P.C );                                               	title('Covariance');            axis square
% subplot(2, 2, 3); plot( t, P.V(:,1), 'k', t, P.V(:,2), 'b', t, P.V(:,3), 'r');	title('Principal Components');  axis square   
% subplot(2, 2, 4); plot3( P.Proj(:, 1), P.Proj(:,2), P.Proj(:,3), '.' );        	title('Projection');            axis square; grid on; rotate3d on



% get covariance matrix of the TRANSPOSE of spike array (waveforms need
% to be in the rows for cov to give what we want)
P.C = cov( data ); 

% get eigenvectors & eigenvalues - these are pre-sorted in order of
% ASCENDING eigenvalue
[P.V, P.D] = eig( P.C );
eigvals = diag(P.D);

% sort in order of DESCENDING eigenvalues
[eigvals, indx] = sort( eigvals , 'descend');
P.V = P.V(:, indx);

% Project original waveforms into eigenvector space
P.Proj = data * [P.V];
    
    
    
    
    
    
    
    
    
