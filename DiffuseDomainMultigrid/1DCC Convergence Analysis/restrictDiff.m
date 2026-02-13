function cDeW = restrictDiff(DeW)

nfEw = length(DeW);
nf = nfEw-1;
nc = nf/2;
cDeW = zeros(nc+1,1);

cDeW(1:nc+1) = DeW(1:2:nf+1);

end