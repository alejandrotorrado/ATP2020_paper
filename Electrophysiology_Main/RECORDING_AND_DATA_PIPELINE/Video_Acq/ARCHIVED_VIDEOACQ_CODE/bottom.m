function bottom(obj, event, string_arg, save_int) % obj, event


global handles
disp(['Count: ' num2str(handles.count)]);
handles.count = handles.count+1;

handles.frametimes(end+1,1) = unixtime(clock);

tic
frame = snapshot(handles.cam);
frame = rgb2gray(frame);
frame = imresize(frame,0.5);


writeVideo(handles.v,frame);

toc

% write new lines to the mac mini stat file on the raspi - it will get
% angry if this stops increasing in time. 
% if unixtime(clock) - handles.timezero>50; % only do this after ~50s of data
%     uploadchecker_pi(unixtime(clock),handles.checkerDir,0);
%     handles.timezero = unixtime(clock);
% end
% % - - - - - - - - - - - - - -


% if vid is 15min, rename prior file to match (parsing)
if ~mod(handles.count,save_int);
    close(handles.v)
    
    ggg = imagesc(frame);
    handles.fileCount = handles.fileCount + 1;
    
    fname = ['Frametimes_' num2str(handles.fileCount) '.mat'];
    frametimes = handles.frametimes;
    save([handles.vDir  fname],'frametimes');
    
    handles.frametimes = [];
    
    term_cmd = ['mv ' handles.vDir filesep handles.vFile ' ' handles.vDir filesep handles.vFile(1:end-4) '_' num2str(handles.fileCount) '.avi' ];
    
    changed_name = unix(term_cmd);
    open(handles.v);
    
end


