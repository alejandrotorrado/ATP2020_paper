function [MSE,normMSE,RMSE,meanpeak] = get_WF_MSE(wf1,wf2,plot_verbose)
% get_WF_MSE
%
% ATP, 1/18/2018
%
% This function compares two waveforms and returns the MSE (mean squared
% error), RMSE (root-mean-square error, sqrt(MSE)), and peak-normalized
% MSE. The two waveforms will usually be average spike waveforms.
%
% INPUTS
%
% __wf1 : first waveform, nx1 or 1xn array.
% __wf2 : second waveform, nx1 or 1xn array.

% optional argument
if nargin < 3
    plot_verbose = 0;
end

% check for input errors
if size(wf1,1) == 1
    wf1 = wf1';
end
if size(wf2,1) == 1
    wf2 = wf2';
end
if size(wf1,1) ~= size(wf2,1)
    error('get_WF_MSE: the two waveforms should have equal dimensions.\n');
end
if all(size(wf1) > 1) || all(size(wf2) > 1)
    error('Input waveforms should be 1xN or Nx1 arrays.');
end

% peak-align
[swa,ewa,D] = alignsignals(wf1,wf2);

% cut
if size(swa,1) > size(ewa,1)
    sw = swa(1+abs(D):end-abs(D));
    ew = ewa(1+abs(D):end);
elseif size(ewa,1) > size(swa,1)
    sw = swa(1+abs(D):end);
    ew = ewa(1+abs(D):end-abs(D));
elseif size(ewa,1) == size(swa,1)
    sw = swa;
    ew = ewa;
end

if plot_verbose
    figure(123182);
    set(gcf,'NumberTitle','off');
    plot(sw,'linewidth',2);
    hold on;
    plot(ew,'linewidth',2);
    title('Close figure to continue.','fontsize',16);
    uiwait(gcf);
    
end

% compare
meanpeak = mean([abs(min(sw)),abs(min(ew))]);
MSE = mean((sw - ew).^2);
RMSE = sqrt(MSE);
normMSE = MSE / meanpeak;
