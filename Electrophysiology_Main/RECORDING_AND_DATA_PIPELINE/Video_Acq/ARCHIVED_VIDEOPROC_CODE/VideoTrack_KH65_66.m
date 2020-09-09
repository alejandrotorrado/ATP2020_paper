% movement tracking main code

clear all, close all, clc
recname = input('Recording name (e.g. KH67_68) ?  ','s');
anim1 = input('Animal 1 name (e.g. KH67) ?  ','s');
anim2 = input('Animal 2 name (e.g. KH68) ?  ','s');
acqdir = fullfile('F:\Animal VIDEO frames',recname);

if ~exist(acqdir)
    err_str = sprintf('VIDEOTRACK : could not find directory %s. Make sure recording name (%s) is correct.\n',...
        acqdir,recname);
    error(err_str);
end

savedir = fullfile(acqdir,'STITCH_VID');
outputdir = fullfile(acqdir,'TRACKING');

if ~exist(savedir,'dir')
    mkdir(savedir);
end

avifiles = dir([acqdir filesep '*.avi']);
matfiles = dir([acqdir filesep '*.mat']);

n_avi = size(avifiles,1);
n_csv = size(matfiles,1);

% index AVI files
avicharlist = char(avifiles.name);
avilist = cellstr(avicharlist);
ext = cell2mat(regexp(avilist,'.avi'));
underscore = cell2mat(regexp(avilist,[recname '_'],'end'));
for cc = 1:size(avicharlist,1)
    temp = cellstr(avicharlist(cc,underscore(cc,1)+1:ext(cc,1)-1));
    avi_index(cc,1) = cc;
    avi_index(cc,2) = str2double(temp{1});
end

avi_index = sortrows(avi_index,2);

% index CSV files
matcharlist = char(matfiles.name);

matlist = cellstr(matcharlist);
ext = cell2mat(regexp(matlist,'.mat'));
underscore = cell2mat(regexp(matlist,['Frametimes_'],'end'));
for cc = 1:size(matcharlist,1)
    temp = cellstr(matcharlist(cc,underscore(cc,1)+1:ext(cc,1)-1));
    mat_index(cc,1) = cc;
    mat_index(cc,2) = str2double(temp{1});
end

mat_index = sortrows(mat_index,2);

for ff = 182:n_avi
    fprintf('Processing file %u out of %u.\n',ff,n_avi);
    temp_row = find(avi_index(:,2)==ff);
    if ~isempty(temp_row)
        temp_idx = avi_index(temp_row,1);
        new_avi_name = fullfile(acqdir,avifiles(temp_idx).name);
        temp_avi = VideoReader(new_avi_name);
        
        mat_name = fullfile(acqdir,matfiles(temp_idx).name);
        
        if ff == 1 % if first file, get ROIs
            read_new_avi = VideoReader(new_avi_name);
            firstframe = readFrame(read_new_avi);
            [ROI_1, ROI_2, corners_1, corners_2] = selectMasks(firstframe);
        end
        
        % now do the tracking
        movementTracking_byChunk_fxn_wMat(new_avi_name,ff,outputdir,...
            anim1,anim2,mat_name,ROI_1,ROI_2,corners_1,corners_2);
    end
end


fprintf('\n\n*** All done! ***\n\n');

% now stitch together tracking output for each chunk






