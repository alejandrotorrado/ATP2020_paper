function [out_path] = add_anaconda_to_path(anaconda_path)
% Add anaconda python and libraries path to current environment - this is
% required to be able to use python callers


if nargin == 1
    % if path is provided, just add it to matlab path
    fprintf(['\nAnaconda path provided by user:\n%s',...
        'Adding to path.\n\n'],anaconda_path);
    
    new_path = [anaconda_path ';' getenv('PATH')];
    fprintf('New PATH variable: %s\n\n',new_path);
    
else
    % if not, try to find it
    fprintf(['\nAnaconda path not provided by user. Looking for it.\n',...
        'This will take a few minutes. Be patient please...\n\n']);
    
    t0 = tic;
    [~,cmdout] = system('find "/" -name "anaconda" -type d 2>&1 | grep -v "Permission denied"');
    t1 = toc(t0);
    fprintf('That took %.2f seconds.\n\n',t1);
       
    cmdsplit = regexp(cmdout,'\n','split');
    
    cmdanac = cellfun(@(x) regexp(x,'anaconda'),cmdsplit,'uniformoutput',0);
    
    new_path = getenv('PATH');
    for aa = 1:length(cmdanac)
        if ~isempty(cmdanac{aa})
            str_to_add = [cmdsplit{aa} filesep 'bin'];
            fprintf('Found anaconda folder: %s\nAdding to path.\n\n',str_to_add);
            new_path = [str_to_add ':' new_path];
            fprintf('New PATH variable: %s\n\n',new_path);
        end
    end
    
end

fprintf('\nRemember to add newly created PATH string to your current Matlab path!\n\n');

out_path = new_path;