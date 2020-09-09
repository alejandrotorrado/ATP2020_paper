function [time_array, file_list] = getCsvFileTimes(csvdir)

allcsv = dir([csvdir filesep '*.csv']);

t0 = tic;
for f = 1:size(allcsv,1)
    
    dt = allcsv(f).date;
    d = datetime(dt);
    ut = posixtime(d);
    temp_arr(f,1) = ut;
    temp_arr(f,2) = f;
    
end
t00 = toc(t0);
fprintf('Getting csv file times took %.2f seconds.\n',t00);

time_array = sortrows(temp_arr,1);
file_list = allcsv;



