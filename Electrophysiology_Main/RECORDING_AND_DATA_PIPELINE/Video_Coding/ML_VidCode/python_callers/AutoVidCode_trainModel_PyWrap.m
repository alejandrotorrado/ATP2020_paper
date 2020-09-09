function [pyMdl_path,OOBerr,saveMdl_path] = AutoVidCode_trainModel_PyWrap(train_set,anim,mlmode,nTrees,pc_name)
% This is the function that takes the user-coded training data set and
% trains a ML model with it.
% INPUTS:
% - train_set: this is the training data set
% - anim: animal ID name
% - mlmode: string specifying which algorithm to use. Current options:
%           * 'RF' : Random Forest
%           * 'ECOC' : error-correcting output codes
% - nTrees: number of decision trees to use for RF model training
%
% Author: ATP, June 2017

%% SETUP;
% keyboard;
% set up default arguments

if nargin < 4
    nTrees = 200;
end
if isempty(nTrees)
    nTrees = 200;
end

% this is set up to work with the Turrigiano Lab server found at:
% files.brandeis.edu/turrigiano-lab
%
% Should work as long as server is mounted on Mac, on windows be sure to
% mount server on drive Z: (or change path below)
if ismac
    % hard-coded repository of trained models
    base_dir = '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Video_Coding/ML_VidCode/TRAINED_MODELS';
    % hard-coded repository for python scripts called by matlab
    pywrap_dir = '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Video_Coding/ML_VidCode/python_callers';
elseif ispc
    base_dir = 'Z:\RECORDING_AND_DATA_PIPELINE\Video_Coding\ML_VidCode\TRAINED_MODELS';
    pywrap_dir = 'Z:\RECORDING_AND_DATA_PIPELINE\Video_Coding\ML_VidCode\python_callers';
end
% pywrap_dir = '/Volumes/turrigiano-lab/ATP_MAIN/CODE/VidCode_MLS/python_callers';
% pywrap_dir = '/Users/atorrado/Desktop/python_callers';

% make directory for this animal where trained model will be saved
mdl_dir = [base_dir filesep mlmode filesep anim];
if ~exist(mdl_dir,'dir'), mkdir(mdl_dir); end

fprintf('Training %s model.\n',mlmode);

% check which ML algorithm to use
switch mlmode
    
    case 'RF'
        %% RANDOM FOREST
        % first save file to be able to read it into python
        training = train_set(:,1:end-1);
        scoredclass = train_set(:,end);
        homefolder = getuserdir;
%         keyboard;
        if isempty(homefolder)
            homefolder = 'C:\Users\Admin\';
        end
        if strcmp(pc_name,'Jerboa')
            tempdir = [homefolder filesep 'Documents' filesep 'MATLAB' filesep 'tempPyVidcodeFiles'];
        elseif strcmp(pc_name,'Malabar')
            tempdir = [homefolder filesep 'MATLAB' filesep 'tempPyVidcodeFiles'];
        end
        if ~exist(tempdir, 'dir'), mkdir(tempdir); end
        tempfile_code = 'temp_train_file.mat';
        tempfile = fullfile(tempdir,tempfile_code);
        if exist(tempfile,'file'), delete(tempfile); end
        save(tempfile,'training','scoredclass');
        
        rf_t0 = tic;
        
        
        %% here call Python code instead of Matlab treebagger
        % after call to code, load OOB error from OOB_error.mat in TEMP
        % folder
        sys_py_call = ['python ' pywrap_dir filesep ...
            'AutoVidCode_trainModel_Python.py -f ' tempfile ' -nt ' num2str(nTrees)...
            ' -sd ' tempdir];

        py_train_status = system(sys_py_call);
        OOB_file_code = 'OOB_error.mat';
        try
        oob_file_load = load([tempdir filesep OOB_file_code]);
        catch
            keyboard;
        end
        OOB_Error = oob_file_load.OOB_error;
        
        
        % this is old
        %{
        Mdl = TreeBagger(nTrees,training,scoredclass,'OOBPrediction','On',...
            'Method','classification', 'PredictorNames',...
            {'D','T','diff','mvtZ','emgZ','emgV','slDT','slD','slT','slMv','prev'},...
            'OOBPredictorImportance','On'); %'OOBPredictorImportance','On'
        OOB_Error_Ensemble = oobError(Mdl);
        OOB_Error = OOB_Error_Ensemble(end);
        %}
        
        rf_t1 = toc(rf_t0);
        fprintf('Time elapsed: %.2f.\n\n',rf_t1);
        
    case 'ECOC'
        %% ERROR CORRECTING OUTPUT CODE
        fprinf('ECOC not available yet. Sorry! Try RF.\n');
        
end

pyMdl_path = [tempdir filesep 'pyMdl.pkl'];
OOBerr = OOB_Error;
saveMdl_path = [mdl_dir filesep mlmode '_pyMdl.mat'];
        