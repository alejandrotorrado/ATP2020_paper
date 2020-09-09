
function makeSpikeFiles(blockpath,animalName,filedate,chunk,LocDir,ANIM,recstart,skipthesechans,T_CHAN,LGN_EXP)

% Based on analyzemystuff_v10 - ATP, September 2016
% This function uses the TDT raw data (SEV format) and creates both the
% .mat spike files (1 per block) and the binary spike files (1 per channel)
% that will be used in later analysis.
%
% INPUT ARGUMENTS
% __blockpath: directory path for the current block. (string)
% __animalName: animal identifier (string)
% __filedate: date for experiment start day, format mmddyy (string)
% __chunk: current block number (integer)
% __LocDir: local directory on which to store spike files before backup
%           (string)
% __ANIM: numerical animal identifier, i.e. animal 1 or 2 (integer)
% __recstart: TSQ timestamp for start of current recording chunk (double,
%             unix timestamp)
% __skipthesechans: array of channels to skip during spike extraction.
%                   Usually these are EMG and reference channels (array of 
%                   doubles)
% keyboard;

% check if T_CHAN argument has been passed, otherwise set to default
switch nargin 
    case 8
        T_CHAN = 66;
        LGN_EXP = 0;
    case 9
        LGN_EXP = 0;
end


% allspikes directory
spikedir = [LocDir filesep animalName filesep 'ALLSPIKES'];

% create log filename (will pass to spikesFromTDTRaw_v10)
spikeLog = [spikedir filesep animalName '_SpikeData_LOG.txt'];

% Get start time of block from dedicated timer channel
time0 = SEV2mat(blockpath,'channel',T_CHAN,'eventname','RAWX','VERBOSE',0);
% make sure to convert to double! If you don't, you will get timing errors
% that are a VERY ANNOYING to deal with.
try
    dat1 = double(time0.RAWX.data(1));
catch
    keyboard;
end
gotime = dat1 + recstart;


% set up inputs to binary spike/time creation function
chans = 1:32;
goodchans = setdiff(chans,skipthesechans);


t01 = tic;
% get spikes from raw data (filtering, thresholding, etc... is done here)
[timesA, spikesA, INFOA, timesB, spikesB, INFOB] = spikesFromTDTRaw_Upsample(blockpath,chunk,ANIM,gotime,spikeLog,goodchans,LGN_EXP);
t11 = toc(t01);
fprintf('%s: Took %.2f seconds to extract, align and interpolate spikes.\n',animalName,t11);




if ANIM==1 % Save the data from animal A, ch 1-32
    disp('Saving spike.mat file.');
    save([spikedir filesep (animalName) '-' (filedate) '-' num2str(chunk) '-' 'spikes.mat'],...
        'spikesA','timesA','INFOA','-v7.3');
    
    if all(cellfun(@isempty,spikesA))
        disp('No spikes found! Skipping binary file update.');
    else
        nSAMPLES = 0; % initialize as 0
        nsamp_counter = 0;
        while nSAMPLES == 0 % while still 0, cycle thru channels to get it
            nsamp_counter = nsamp_counter + 1;
            nSAMPLES = size(spikesA{nsamp_counter},2);
        end
        
        % Create binary spike files for this animal
        try
            t02 = tic;
            makeSpikeBinaries(spikesA,timesA,INFOA,spikedir,chans,goodchans,animalName,nSAMPLES);
            t12 = toc(t02);
            fprintf('%s: Took %.2f seconds to append spikes to binary file.\n',animalName,t12);
        catch
            warning('makeSpikeFiles: error when trying to run makeSpikeBinaries. Keyboarded.');
            keyboard;
        end
        disp('Done writing binary data.');
    end
    
end

if ANIM==2 % Save the data from animal B, ch 33-64
    
    
    disp('Saving spike.mat file.');
    save([spikedir filesep (animalName) '-' (filedate) '-' num2str(chunk) '-' 'spikes.mat'],...
        'spikesB','timesB','INFOB','-v7.3');
    
    if all(cellfun(@isempty,spikesB))
        disp('No spikes found! Skipping binary file update.');
    else
        nSAMPLES = 0; % initialize as 0
        nsamp_counter = 0;
        while nSAMPLES == 0 % while still 0, cycle thru channels to get it
            nsamp_counter = nsamp_counter + 1;
            nSAMPLES = size(spikesB{nsamp_counter},2);
        end
        
        % Create binary spike files for this animal
        try
            t02 = tic;
            makeSpikeBinaries(spikesB,timesB,INFOB,spikedir,chans,goodchans,animalName,nSAMPLES);
            t12 = toc(t02);
            fprintf('%s: Took %.2f seconds to append spikes to binary file.\n',animalName,t12);
        catch
            warning('makeSpikeFiles: error when trying to run makeSpikeBinaries. Keyboarded.');
            keyboard;
        end
        disp('Done writing binary data.');
    end
end


end



