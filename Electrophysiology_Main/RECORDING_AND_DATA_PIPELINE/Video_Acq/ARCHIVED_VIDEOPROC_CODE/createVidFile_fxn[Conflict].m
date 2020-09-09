function [vidfile,vidfilename,vidTimestamps,frametimes] = ...
    createVidFile_fxn(rawVideoDir,saveDir,animal,framerate)

% createVidFile_fxn
%
% Function version of the createVidFile script
%
% This takes the raw video files saved by the MacMini and creates a single
% avi file containing the whole video recording. This will then be used to
% perform motion tracking using blob analysis.
%
% INPUTS
%   rawVideoDir: directory containing the raw video files created by the
%                MacMini.
%   saveDir: directory in which to save the processed video (.avi file) and
%            the movement data that will be calculated later.
%   animal: animal name (e.g. KH50 or KH60_61). This will be prepended to
%           the filename of the .avi video and movement file so should be
%           a unique identifier for that recording.
%   framerate: arbitrary framerate for the .avi video file. Not important
%              but required. Generally set to 30.
%
% OUTPUTS
%   vidfile: videoWriter object pointing to the .avi video file created by 
%            stitching together all of the frames recorded by the MacMini.
%   vidfilename: filename of the .avi video file.
%   vidTimestamps: filename of the .mat file containing the timestamps of
%                  each frame in the .avi video file
%   frametimes: array containing unixtime timestamps for each frame in the
%               .avi video file
%
%
%
% ATP, April 2016

% use these for testing:
% rawVideoDir = uigetdir(cd,'raw video');
% saveDir = uigetdir(cd,'save dir');
% animal = input('animal:  ','s');
% framerate = 30;

vidTimestamps = [saveDir filesep animal '_frameTimes.mat'];
vidfilename = [saveDir filesep animal '_fullVideo.avi'];
vidfile = VideoWriter(vidfilename);
vidfile.FrameRate = framerate;


rawfiles = dir([rawVideoDir filesep '*.mat']);
N_files  = size(rawfiles,1);
filedates(:,1) = datenum(cat(1,rawfiles.date));
filedates(:,2) = 1:N_files;
filedates = sortrows(filedates,1);


open(vidfile);
t00 = tic;

framecounter = 0;
for ii = 1:N_files
    t10 = tic;
    imgindex = filedates(ii,2);
    imgfile = [rawVideoDir filesep rawfiles(imgindex).name];
    tic
    dataloader = load(imgfile,'saveframe');
    toc
    images = dataloader.saveframe.image;
    times  = dataloader.saveframe.time;
    
    for ee = 1:size(images,2)
        framecounter = framecounter + 1;
        
        img = images{ee};
        
        if ~isempty(img)
            frametimes(framecounter,1) = times(ee);
            writeVideo(vidfile,img);
        else
            fprintf('Skipping empty frame!');
            framecounter = framecounter - 1;
        end
        
    end
    t11 = toc(t10);
    fprintf('  Took %3.2f seconds to write file %u.\n\n',t11,ii);
end
close(vidfile);
t01 = toc(t00);

fprintf('\n Took %3.2f seconds to write full video (%u files, ~%u hours).\n',t01,N_files,N_files/4);

fprintf('Now saving video frame timestamps to %s...\n',vidTimestamps);
tic
save(vidTimestamps,'frametimes');
toc
fprintf('Done!\n\n');
