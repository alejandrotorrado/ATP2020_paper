function [out_Mdl,out_err,Mdl_Path] = AutoVidCode_trainModel(train_set,anim,mlmode,nTrees,display_flag)
% This is the function that takes the user-coded training data set and
% trains a ML model with it.
% INPUTS:
% - train_set: this is the training data set
% - anim: animal ID name
% - mlmode: string specifying which algorithm to use. Current options:
%           * 'RF' : Random Forest
%           * 'ECOC' : error-correcting output codes
% - nTrees: number of decision trees to use for RF model training
% - display_flag: 1 to display error predictions, 0 to not display (default)
%
% Author: ATP, June 2017

%% SETUP;

% set up default arguments
if nargin < 4
    nTrees = 200;
end
if nargin < 5
    display_flag = 0;
end
if isempty(nTrees)
    nTrees = 200;
end

% hard-coded repository
base_dir = '/Volumes/turrigiano-lab/RECORDING_AND_DATA_PIPELINE/Video_Coding/ML_VidCode/TRAINED_MODELS';

% make directory for this animal where trained model will be saved
mdl_dir = [base_dir filesep mlmode filesep anim];
if ~exist(mdl_dir,'dir'), mkdir(mdl_dir); end

fprintf('Training %s model.\n',mlmode);

% check which ML algorithm to sue
switch mlmode
    
    case 'RF'
        %% RANDOM FOREST
        training = train_set(:,1:end-1);
        scoredclass = train_set(:,end);
        
        rf_t0 = tic;
        
        Mdl = TreeBagger(nTrees,training,scoredclass,'OOBPrediction','On',...
            'Method','classification', 'PredictorNames',...
            {'D','T','diff','mvtZ','emgZ','emgV','slDT','slD','slT','slMv','prev'},...
            'OOBPredictorImportance','On'); %'OOBPredictorImportance','On'
        OOB_Error_Ensemble = oobError(Mdl);
        OOB_Error = OOB_Error_Ensemble(end);
        
        rf_t1 = toc(rf_t0);
        fprintf('Time elapsed: %.2f.\n\n',rf_t1);
        
        if display_flag
            
            oobfig = figure(); hold on;
            title('Close this figure to continue running code');
            set(oobfig,'units','normalized','position',[.1 .1 .6 .5]);
            oobErrorBaggedEnsemble = oobError(Mdl);
            plot(oobErrorBaggedEnsemble,'color',[0 0 0 .7],'linewidth',2)
            xlabel 'Number of grown trees';
            ylabel 'Out-of-bag classification error';
            set(gca,'ylim',[0 0.4],'fontsize',12,'Xcolor','k','YColor','k');
            [minoob,minoob_idx] = min(oobErrorBaggedEnsemble);
            minline = refline(0,minoob);
            minline.Color = [0.30 .75 .93 0.7];
            minline.LineWidth = 2;
            text(40,0.35,sprintf('Best performance (OOB Error = %.2f%%) with %u trees.',100*minoob,minoob_idx),...
                'fontsize',11);
            text(40,0.31,sprintf('Final performance: OOB Error = %.2f%%',100*oobErrorBaggedEnsemble(end)),...
                'fontsize',11);
            
            uiwait(oobfig);
            
        end
        
    case 'ECOC'
        %% ERROR CORRECTING OUTPUT CODE
        fprinf('ECOC not available yet. Sorry! Try RF.\n');
        
end

Mdl_Path = [mdl_dir filesep mlmode '_Mdl.mat'];
out_Mdl = Mdl;
out_err = OOB_Error;
        