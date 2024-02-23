function hz = forceSquarePlateRandomDistributed(x,y,z)
L1 = 1; %pi;
L2 = 1; %exp(1);

d1 = L1/20;
d2 = L2/20;

hz = zeros(size(x));
idxx = find(x<= 9*L1/10 & x>=L1/10 -d1);
idxy = find(y<= L2/2 +d2 & y>=L2/2 -d2);
idxz = find(z);
idx = intersect(idxx, idxy);
idx = intersect(idx,idxz);
hz(idx) = -1 +2*rand(size(x));
end