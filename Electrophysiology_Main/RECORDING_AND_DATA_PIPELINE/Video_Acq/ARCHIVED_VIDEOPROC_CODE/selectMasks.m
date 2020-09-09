function [roi1,roi2,region_1,region_2] = selectMasks(vidframe)

imshow(vidframe);
set(gcf,'name','Select corners of first ROI (animal 1). Close figure to continue.',...
    'numbertitle','off');
msg1 = msgbox({'Select corners of first ROI (top left, then bottom right).' ...
    'Close this figure when you are finished to continue to the second ROI selection.'},...
    'ATTENTION!','modal');

for cc=1:2
    [x(cc,1),y(cc,1)]=ginput(1);
    x(x<1)=1;
    x(x>size(vidframe,2)) = size(vidframe,2);
    y(y<1)=1;
    y(y>size(vidframe,1)) = size(vidframe,1);
end
rectangle('Position',[x(1) y(1) x(2)-x(1) y(2)-y(1)],'EdgeColor','r','linewidth',1.5);
region_1 = round([x y]);

% make binary image mask
roi1 = zeros(size(vidframe,1),size(vidframe,2));
roi1(region_1(1,2):region_1(2,2),region_1(1,1):region_1(2,1)) = 1;


uiwait(gcf);

imshow(vidframe);
set(gcf,'name','Select corners of second ROI (animal 2). Close figure to continue.',...
    'numbertitle','off');
msg2 = msgbox({'Select corners of second ROI (top left, then bottom right).' ...
    'Close this figure when you are finished to continue to the second ROI selection.'},...
    'ATTENTION!','modal');

for cc=1:2
    [x(cc,1),y(cc,1)]=ginput(1);
    x(x<1)=1;
    x(x>size(vidframe,2)) = size(vidframe,2);
    y(y<1)=1;
    y(y>size(vidframe,1)) = size(vidframe,1);
end
rectangle('Position',[x(1) y(1) x(2)-x(1) y(2)-y(1)],'EdgeColor','r','linewidth',1.5);
region_2 = round([x y]);

% make binary image mask
roi2 = zeros(size(vidframe,1),size(vidframe,2));
roi2(region_2(1,2):region_2(2,2),region_2(1,1):region_2(2,1)) = 1;

uiwait(gcf);

imshow(vidframe);
msg3 = msgbox({'These are the ROIs you selected for animals 1 (yellow) and animal 2 (blue)' ...
    'Close this figure to continue.'},...
    'FINAL CHECK!','modal');
rectangle('Position',[region_1(1,1) region_1(1,2) region_1(2,1)-region_1(1,1) region_1(2,2)-region_1(1,2)],'EdgeColor',[0.98 0.64 0.10],'linewidth',1.5);
rectangle('Position',[region_2(1,1) region_2(1,2) region_2(2,1)-region_2(1,1) region_2(2,2)-region_2(1,2)],'EdgeColor',[0.30 0.75 0.93],'linewidth',1.5);
uiwait(gcf);



