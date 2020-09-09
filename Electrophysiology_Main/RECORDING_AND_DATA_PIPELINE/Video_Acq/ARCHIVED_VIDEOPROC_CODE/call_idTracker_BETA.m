% RUN idTracker
%
% ATP, April 2016
% Modifying idTracker.m
% Source code by:
% Pérez-Escudero, Alfonso, et al. "idTracker: tracking individuals in a
% group by automatic identification of unmarked animals." Nature methods
% 11.7 (2014): 743-748.
%
% This is the function you need to run:
% idTracker_ATP(directorio,nombrearchivo,directorio_destino,n_peces,umbral,reutiliza,roi,cambiacontraste,referencias,mascara_intensmed,solodatosegm,nRefFrames,minBlobPixels,reduceResolution)
%
% Argument explanation:
% __directorio: video file directory
% __nombrearchivo: video file name
% __diretorio_destino: saving directory
% __n_peces: number of animals to track
% __umbral: intensity threshold for blob detection
% __reutiliza: don't care about this
% __roi: ROI for tracking. Will get this from first frame
% __cambiacontraste: don't care about this
% __referencias: don't care about this
% __mascara_intensmed: don't care about this
% __solodatosegm: don't care about this
% __nRefFrames: number of reference frames (for bkg & id tracking)
% __minBlobPixels: minimum area of blobs to detect
% __reduceResolution: reduce video res (for speed)
%
% All arguments we don't care about will be set to their default values as
% given in the source code. All others will be set before calling the
% function.

clear all, close all, clc

%% CHOOSE VIDEO FILE AND SAVE DIRECTORY
[video_file, video_dir] = uigetfile('*.*', 'Pick your first video file.');
save_dir = [video_dir 'segm' filesep];
if ~exist(save_dir,'dir'), mkdir(save_dir); end


%% SELECT ROI FROM FIRST FRAME
vidSrc = VideoReader([video_dir  video_file]);

% get number of frames in vid
n_frames = get(vidSrc,'numberOfFrames');

firstframe = read(vidSrc,1);
firsftrame = rgb2gray(firstframe);

imshow(firstframe);
set(gcf,'name','Select corners of ROI. Close figure to continue.',...
    'numbertitle','off');
msg1 = msgbox({'Select corners of ROI (click top left, then bottom right).' ...
    'Close this figure when you are finished to continue.'},...
    'ATTENTION!','modal');

for cc=1:2
    [x(cc,1),y(cc,1)]=ginput(1);
    x(x<1)=1;
    x(x>size(firstframe,2)) = size(firstframe,2);
    y(y<1)=1;
    y(y>size(firstframe,1)) = size(firstframe,1);
end 
rectangle('Position',[x(1) y(1) x(2)-x(1) y(2)-y(1)],'EdgeColor','r','linewidth',1.5);
ROI = round([x y]);
pause(1);
close(gcf);

%% GET MEAN INTENSITY FOR ALL FRAMES AND SET INTENSITY THRESHOLD FOR L/D
fprintf('\n\nReading frames for mean intensity calculation.\n');
t0 = tic;
n_chunks = 12;
chunk_sz = 1e3;
chunk_jump = 1e5;
for ee = 1:n_chunks
    chunk{ee} = read(vidSrc,chunk_jump*(ee-1)+[1 chunk_sz]);
    fprintf('Done with %u out of %u chunks...\n',ee,n_chunks);
end
t1 = toc(t0);
fprintf(['\nIt took %.1f seconds to read %u frames, in chunks of %u frames ',...
'spaced %u frames apart.\n\n'],t1,n_chunks*chunk_sz,chunk_sz,chunk_jump);

allframes = cat(4,chunk{:});

mask = zeros(size(allframes,1),size(allframes,2));
mask(ROI(1,2):ROI(2,2),ROI(1,1):ROI(2,1)) = 1;

% avgIntensity = squeeze(mean(mean(mean(allframes(a),1),2),3));
t2 = tic;
for xx = 1:size(allframes,4)
%     disp(xx)
    temp_frame = rgb2gray(allframes(:,:,:,xx));
    avgIntensity(xx) = mean(mean(temp_frame));
    clear temp_frame
end
t3 = toc(t2);
fprintf('Took %.1f seconds to calculate mean intensity for %u frames.\n',t3,n_chunks*chunk_sz);
int_fig = figure();
plot(avgIntensity);
set(gcf,'numbertitle','off');
title('Select intensity threshold for L/D discrimination.');
[~,thresh] = ginput(1);
pause(0.5);
close(int_fig);
fprintf('Intensity threshold set at %.2f.\n\n',thresh);



%% SET OTHER PARAMETERS


% corresponding argument given in comments
n_animals = 1;                          % n_peces
int_thresh_light = 0.70;                % umbral for light frames
int_thresh_dark = 0.85;                 % umbral for dark frames
reuse = false;                          % reutiliza
change_contrast = false;                % cambiacontraste
references = [];                        % referencias
mask_intens = [];                       % mascara_intensmed
solodatosegm = false;                   % solodatosegm
n_ref_frames = n_frames;                % nRefFrames - I set this to n_frames because it doesn't seem to affect bkg calculation 4/15/16
min_blob_area_light = 500;              % min area of blobs for light frames
min_blob_area_dark = 2000;              % min area of blobs for dark frames
max_blob_area_light = 4000;             % max area of blobs for light frames
max_blob_area_dark = 15000;             % max area of blobs for dark frames
intens_thresh = thresh;                 % intensity threshold
res_reduction = 4;                      % reduceResolution
n_bkg_frames = round(n_frames*.1,-5);   % number of frames to use for bkg calculation
disp(n_bkg_frames);
max_nbkg = 400000;
if n_bkg_frames > max_nbkg
    n_bkg_frames = max_nbkg;            % set max number of frames for bkg calc
end

idTracker_ATP(video_dir, video_file, save_dir, n_animals, int_thresh_light,...
    int_thresh_dark, reuse, ROI, change_contrast, references, mask_intens,...
    solodatosegm,n_ref_frames, min_blob_area_light, min_blob_area_dark,...
    max_blob_area_light, max_blob_area_dark, intens_thresh, res_reduction,...
    n_bkg_frames);



