function [folder_cells,folder_list,n_folders] = get_folder_info(topDir,selectName)

% May 2016
% Updated this script for new RS4 naming convention.

direc_ix = dir(topDir);
for xx = 1:length(direc_ix)
    folder_cells{xx,1} = direc_ix(xx).name; % get all folder names
end
% remove all folder names without our "selectName"
dkill   = regexp(folder_cells,selectName,'once');
dkillem = find(cellfun(@isempty,dkill));
folder_cells(dkillem) = [];
tankstr = length(selectName);
% convert to a matrix of folder names for easy access
try
    folder_list = cell2mat(folder_cells);
catch
    disp(['could not convert dfolder to matrix']);
    keyboard;
end
% sort folders by creation time
folder_list = sortrows(folder_list,[tankstr+4,tankstr+5,tankstr+6,tankstr+7,...
    tankstr+9,tankstr+10,tankstr+11,tankstr+12,tankstr+13,tankstr+14]);
% get number of data folders in tank
n_folders = size(folder_cells,1);