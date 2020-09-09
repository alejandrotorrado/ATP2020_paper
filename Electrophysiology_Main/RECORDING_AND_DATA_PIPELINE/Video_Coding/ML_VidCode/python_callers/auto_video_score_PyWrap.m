function [out_score,raw_score,feature_mtx] = auto_video_score_PyWrap(input_data, ml_model_path, tx, bin_size, Fs)
% This function takes in input behavioral data and uses a Machine Learning
% algorithm to automatically score the behavioral state of the animal,
% classifying it as one of 4 states: REM (1), NREM (2), active wake (4),
% quiet wake (5).
%
% 
% Author: ATP
% Date: June 2017
%
% JUNE 2020 NOTE
% A key improvement to the videocoding classifier would be to include more
% past bins as features. As of now the code includes the previous 10-sec
% bin as a feature for the current bin classification. I suspect including
% 2-3 extra bins (i.e. the past 30-40 seconds in total) may improve the
% coding substantially.


% keyboard;

% scoring in 10-sec bins
if nargin < 4
    bin_size = 10; % seconds
end
if nargin < 5
    Fs = 1; % Hz
end

%% COMPUTE BIN EDGES

bin_edges = 0 : bin_size*Fs : size(input_data,1);
if bin_edges(end)~=(size(input_data,1))
    bin_edges(end + 1) = size(input_data,1);
end

% initialize output array
autoscore = zeros(size(bin_edges,2)-1,4);

% find python_callers directory
if ismac
    pywrap_dir = '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Video_Coding/ML_VidCode/python_callers';
elseif ispc
    pywrap_dir = 'Z:\RECORDING_AND_DATA_PIPELINE\Video_Coding\ML_VidCode\python_callers';
end

% state codes
class_labels = [1,2,4,5];

%% LOOP THROUGH BINS
for jj = 1:size(bin_edges,2)-1
    
    % find t0 and t1
    t0      = bin_edges(jj)+1;
    t1      = bin_edges(jj+1);
    autoscore(jj,2) = tx(t0);
    autoscore(jj,3) = tx(t1);
    
    
    %% EXTRACT FEATURES FROM DATA
    
    % Feature list:
    % 1. delta power
    % 2. theta power
    % 3. (delta - theta) difference
    % 4. movement z-score
    % 5. EMG data
    % 6. EMG variance
    % 7. EMG amplitude
    % 8. (delta-theta) slope
    % 9. delta slope
    % 10. theta slope
    % 11. movement slope
    % 12. previous bin score ?
    
    % First slice input array and retrieve each column
    delta_data      = input_data(t0:t1,1);
    theta_data      = input_data(t0:t1,2);
    DTdiff_data     = input_data(t0:t1,3);
    mvt_data        = input_data(t0:t1,4);
    emg_data        = input_data(t0:t1,5);
    
    % Features 1 - 5 are the average of each column of the input matrix
    deltapower      = nanmean(delta_data);
    thetapower      = nanmean(theta_data);
    DT_diff         = nanmean(DTdiff_data);
    mvt_z           = nanmean(mvt_data);
    emg_z           = nanmean(emg_data);
    
    % Features 6 - 11 are computed here
    emg_var         = var(emg_data);
    emg_ampl        = nanmean(abs(diff(emg_data)));
    DT_slope        = (DTdiff_data(end)-DTdiff_data(1))/numel(DTdiff_data);
    delta_slope     = (delta_data(end)-delta_data(1))/numel(delta_data);
    theta_slope     = (theta_data(end)-theta_data(1))/numel(theta_data);
    mvt_slope       = (mvt_data(end)-mvt_data(1))/numel(mvt_data);
    
    % Feature 11 only exists for jj > 1
    if jj > 1
        previous_score = autoscore(jj-1,1);
    else
        % otherwise set to NaN. Make sure algorithm can deal with NaNs!
        % edit: set to 0 instead
        previous_score = 0;
    end
    
    % collect features in array
    feature_array(jj,1)     = deltapower;
    feature_array(jj,2)     = thetapower;
    feature_array(jj,3)     = DT_diff;
    feature_array(jj,4)     = mvt_z;
    feature_array(jj,5)     = emg_z;
    feature_array(jj,6)     = emg_var;
    feature_array(jj,7)     = emg_ampl;
    feature_array(jj,8)     = DT_slope;
    feature_array(jj,9)     = delta_slope;
    feature_array(jj,10)     = theta_slope;
    feature_array(jj,11)    = mvt_slope;
    feature_array(jj,12)    = previous_score;
    
   
end

feature_array(isnan(feature_array)) = 0;

%% FITTING TO MODEL
% save to temp file
homefolder = getenv('HOME');
tempdir = [homefolder filesep 'Documents' filesep 'MATLAB' filesep 'tempPyVidcodeFiles'];
if ~exist(tempdir, 'dir'), mkdir(tempdir); end
tempfeat_code = 'temp_feature_input.mat';
tempfeat = fullfile(tempdir,tempfeat_code);
if exist(tempfeat,'file'), delete(tempfeat); end
save(tempfeat,'feature_array');

% call to python script for scoring
tel0 = tic;
py_score_call = ['python ' pywrap_dir filesep ...
            'auto_video_score_Python.py -f ' tempfeat ' -sd ' tempdir ' -mdl ' ml_model_path];
py_score_status = system(py_score_call);
output_code = 'pyOut.mat';

% now get labels and confidence scores from pyOut.mat in TEMP folder
output_load = load([tempdir filesep output_code]);
labels = output_load.labels';
confidence = output_load.confidence';
label_idxs = (class_labels==labels)';
label_confidence = confidence(label_idxs);

autoscore(:,1) = labels;
autoscore(:,4) = label_confidence;
tel1 = toc(tel0);
fprintf('Done! Time elapsed: %.2f seconds.\n\n',tel1);


%% clear small states
autoscore_out = autoscore;
ticker = 1;
while ticker == 1
    % eliminate states that only last for one 10-sec bin
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % copy the state assignment
    temp = autoscore_out(:,1);
    %find the differential. zero when a state is constant
    temp2 = diff(temp);
    % make it binary. all changes are 1s and constants are 0s
    temp2(temp2~=0) = 1;
    % index the changes
    temp3 = find(temp2 == 1);
    % find the amount of bins between changes
    temp4 = diff(temp3);
    % look for bins changes that are only seperated by 1
    temp5 = find(temp4==1);
    % if you have staggered 1-bin states flanking one another, this will
    % need to run a couple of times to iteratively clean them out. move on
    % when finished:
    if isempty(temp5)
        ticker = 0;
    end
    % get the indices of the bins that have a 1 bin state
    kills = temp3(temp5);
    % overwrite the 1 bin state with the preceding state
    temp(kills+1) = temp(kills);
    % reassign temp into the autocat variable
    autoscore_out(:,1) = temp;
end

% catch bins algorithm didn't score (shouldn't happen with ML?)
if any(unique(autoscore_out(:,1)) == 0)
    disp('Algorithm failed to ID and catch a time bin!');
    keyboard
    % find the zero entry in the first column of autocat, set ee to the row
    % number and then go through the above steps to figure out where it
    % should be ID'ing a state and why it's missing this entry. It might
    % help to turn on "verbose" and see the data to more easily understand
    % what it "should" be identifying.
end


%% output arguments
out_score = autoscore_out;
raw_score = autoscore;
feature_mtx = feature_array;
