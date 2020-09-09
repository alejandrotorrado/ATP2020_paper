rawVideoDir = uigetdir(cd,'Pick the folder with all the raw video files.');
animal = input('Animal(s) name?  ','s');
clearvars -except *Dir animal vidfile framerate
rawfiles = dir([rawVideoDir filesep animal '*.avi']);

for i = 1:size(rawfiles,1)
    name = rawfiles(i).name
    if ~contains(name,num2str(i))==0
        continue;
    else
        keyboard;
    end
end
