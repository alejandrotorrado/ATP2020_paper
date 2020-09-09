function [sum_out] = sum_isi(isidist,bin0,bin1,step)


sum_out = [];
outcount = 0;
for ii = bin0:step:bin1
    outcount = outcount + 1;
    lowbin = ii - step;
    hibin = ii;
    tempisi = sum(isidist<=hibin) - sum(isidist<=lowbin);

    sum_out(outcount) = tempisi;

end