close all; clear all; clc;

TANK = 'EXAMPLE';
BLOCK = 'Block-1';
REF_EPOC = 'Levl';
SNIP_STORE = 'Spik';
SORTID = 'TankSort';
CHANNEL = 1;

data = TDT2mat(TANK, BLOCK, 'TYPE', [2, 3, 5], 'SORTNAME', SORTID);
data = TDTfilter(data, REF_EPOC, 'TIME', [-0.02, 0.07], 'TIMEREF', 1);

i = find(data.snips.Spik.chan == CHANNEL);
TS = data.snips.(SNIP_STORE).ts(i)';

% extract each trial from timestamps
new_trials = find(diff(TS) < 0);
new_trials = [1 new_trials length(TS)];

%plot raster
subplot(2,1,1)
for x = 2:length(new_trials)
    trial = TS(new_trials(x-1)+1:new_trials(x));
    if ~isempty(trial)
        plot(trial, x-1, '.', 'MarkerEdgeColor','k', 'MarkerSize',10)
    end
    hold on;
end
line([0 0], [1, x-1], 'Color','r', 'LineStyle','--')
axis tight;
ylabel('trial number')
xlabel('time, s')
title('Raster')

NBINS = floor(numel(TS)/10);
subplot(2,1,2)
hist(TS, NBINS);
N = hist(TS, NBINS);
hold on;
line([0 0], [0, max(N)*1.1], 'Color','r', 'LineStyle','--')
axis tight;
ylabel('number of occurrences')
xlabel('time, s')
title('Histogram')