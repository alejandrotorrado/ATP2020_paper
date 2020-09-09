numfiles = 4;
fdir = 'D:\AT39_VIDCODE\saving_dir';
filename = 'AT39_STATETIMES_';
close= '.mat';

statetimes = [];

for n=1:numfiles
    a= num2str(n);
    loadvariable = load([fdir filesep filename a close]);
    statetimes = [statetimes; loadvariable.statetimes]; %takes empy cell array, S (defined above), and adds each loadvariable.statetimes file from 1 to 4% 
    
end

save([fdir filesep filename 'FINAL' close], 'statetimes');


