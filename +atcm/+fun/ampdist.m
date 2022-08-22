function y = ampdist(x,N)

[bincounts,binedges] = histcounts(real(x),N);

binedges = binedges(1:end-1)+diff(binedges)/2;

y = [log(bincounts(:)); binedges(:)];

