function [OUTPUT] = alignrows(A,B,xrange,ptscheck)
% ALIGNROWS is a vectorized solution to aligning many trials or samples of
% data (stored in a matrix) on some value that differs by row. For
% example, you may want to extract 10 datapoints from each row, centered
% on the max or min value of that row. ALIGNROWS will return the subset
% matrix that is centered on a vector of nrows indices.
%
% INPUTS:
%   A:
%      A is a matrix in which each row is an event/trial/waveform.
%   B:
%      B is a vector of indices that correspond to a point in each row of
%      A. For example, B might be the location of the minimum value, by
%      row, in A. OUTPUT will be a matrix with each row aligned on the
%      values indicated in B.
%   Range:
%      Range is the range of data to be extracted around the row points
%      indexed in B. Negative values indicate data to the left of B
%      points.
%   Ptscheck:
%       n points -1 in each sample. Assuming that each row is a sample,
%       this is the number of columns minus one (this makes processing
%       spike extraction easier at the expense of making other use a little
%       more confusing). 
% OUTPUT:
%      OUTPUT is a subset of the matrix A, each row centered on the
%      indices written in B.
%

% KBH and ATP 6/23/16

% keyboard;

if isrow(B)
    B = B';
end

if ~isrow(xrange)
    xrange = xrange';
end


[nrow,ncol] = size(A);
if ncol ~= ptscheck+1 && nrow == ptscheck+1
    A = A';
    [nrow,ncol] = size(A);
elseif ~any([ncol nrow] == ptscheck+1)
    disp('STUCK IN ALIGNROWS.M BECAUSE THE DATA INPUT SEEMS WRONG.');
    keyboard
end

% if nrow<ncol;
%     disp('There are more columns than rows in your data. This is unusual.');
%     flp = input('Would you like to rotate the matrix to fix this? Y or N ','s');
%     
%     if strcmp(flp,'Y') || strcmp(flp,'y');
%         A = A';
%     end
% end

xrange       = repmat(xrange,size(B,1),1);
B           = repmat(B,1,size(xrange,2));
C           = B + xrange;

% Get rid of 0 entries
C(C<1)      = 1;
% Get rid of out of range entries
C(C>ncol) = ncol;
% [outrow, outcol] = find(C>ncol);
% C(:,outcol) = [];



% The above line "C(C<1)      = 1;" is the current solution to finding
% multiple column positions of values less than 1.
% D           = C(:,1);
% E           = D(:,1) < 1;
% C(E,1)      = C(E,2);
%
% if sum(C(:) == 0)>0;
%     disp('Caught in alignrows. Figure out how to deal with this instance of zeros');
%     keyboard
% end

clear B D E
% keyboard;
p = 1:size(C,1);
P = repmat(p',1,size(C,2));

clearvars -except A P C xrange
try
    ind = sub2ind(size(A),P,C);
catch
    keyboard
end

clear P C

OUTPUT = A(ind);