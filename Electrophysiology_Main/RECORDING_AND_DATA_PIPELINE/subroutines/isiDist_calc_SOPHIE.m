function [x,y,isiOUT] = isiDist_calc_SOPHIE(stimes,t0,t1,bin)
% 
% This will calculate the percent of spike ISIs less than 2.5 milliseconds 
% in a spike train STIMES between times T0 and T1, operating with a bin 
% window of size BIN.
%
edges = [ t0: bin : t1 ];


Y = discretize(stimes,edges);
b = zeros(2,size(edges,2));


for ee = 1:size(edges,2);
    
    a = diff(stimes(Y == ee));
    b(1,ee) = 100*(sum(a<0.0025)/size(a,1));

end

b(2,:) = edges;

isiOUT = b;
x = b(1,:);
y = b(2,:);
% figure
% plot(b(2,:),b(1,:));