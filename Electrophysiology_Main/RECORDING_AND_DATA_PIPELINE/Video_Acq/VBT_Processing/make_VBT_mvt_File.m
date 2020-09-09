clearIDE

pos_dir = 'E:\JB14_VBT';
csv_dir = 'E:\JB14_VBT';

pos_prefix = 'JB14_';
csv_prefix = 'JB13_14_';

mydir_position = dir([pos_dir filesep pos_prefix 'pos_*.txt']);
num_files_position = length(mydir_position(not([mydir_position.isdir])));
new_file = [];
new_file_ts = [];

%%
bad_counter = 0;
bad_ID = [];

for ii = 1:num_files_position
    fprintf('Processing file %u of %u.\n',ii,num_files_position);
    this_file = [pos_dir filesep mydir_position(ii).name];
    split_1 = regexp(this_file,'_','split');
    file_end = split_1{end};
    split_2 = regexp(file_end,'\.','split');
    file_ID = split_2{1};
    
    csv_file = [csv_dir filesep csv_prefix file_ID '.csv'];
    
    if exist(csv_file,'file')
        % read position file and concatenate in one big matrix
        d1 = dlmread(this_file); 
        % read csv file and concatenate in another matrix
        csv_read = dlmread(csv_file);
        if size(d1,1)== size(csv_read,1)
            if size(d1,2) < 4
                d1 = [d1 repmat(new_file(end,size(d1,2)+1:4),size(d1,1),1)];
            end
           new_file = cat(1, new_file, d1(:,1:4));
           new_file_ts = cat(1, new_file_ts, csv_read);
        else
            fprintf('Found bad file! ID: %s\n',file_ID);
            bad_counter = bad_counter + 1;
            bad_ID{bad_counter} = file_ID;
        end
    else
        fprintf('\tCSV file with ID %s not found. Skipping.\n',file_ID);
    end
end
save('Combined.txt', 'new_file');

d2_a = reducenan(new_file);
d2 = reducezeros(d2_a);
square_root = sqrt(d2(:,3).^2 + d2(:,4).^2);
d3 = [d2 square_root];
POS = d3(:,end);
MVT = abs(diff(POS));
MVT_column = cat(1, NaN , MVT); 
d4 = [d3 POS];
SMOOTH = smooth(MVT_column,20,'moving');

DATA.smooth_movement = SMOOTH';
DATA.frame_times = new_file_ts;
%DATA.mask = 
DATA.nframes = length(d4);

RAW.raw_movement = MVT_column';
RAW.track = d3(:,3:4);

%%
% keyboard;
tic;
disp('Saving data...');
save_file = [fdir filesep pos_prefix 'VBTpymovement.mat'];
save(save_file,'DATA','RAW','-v7.3');
toc;