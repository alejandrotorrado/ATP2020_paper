function [NB] = get_nblocks(rs4,blockdir)

file_index = dir([rs4 filesep blockdir]);
for nn = 1:length(file_index)
    filenames{nn} = file_index(nn).name;
end
% this next part converts to char array and finds maximum 'hour'
% specified in file name, to determine nblocks in curr_dir
filez = char(filenames);
fn_split = regexp(filenames,'_','split');
for ii=1:length(fn_split)
    if length(fn_split{ii}) > 1
        lastone = fn_split{ii}(end);
        dashsplit = regexp(lastone,'-','split');
        if length(dashsplit{1}) == 1
            hourz(ii) = 0;
        else
            dashlast = dashsplit{1}(end);
            digits = regexp(dashlast,'\d*','match');
            hourz(ii) = str2double(digits{1});
        end
    end
end

NB = max(hourz) + 1;

end