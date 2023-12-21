function hz = forceSquarePlateCentered(x,y,z)
L1 = 15; %pi;
L2 = 15; %exp(1);

d1 = L1/20;
d2 = L2/20;

hz = zeros(size(x));
idxx = find(x<= L1/2 +d1 & x>=L1/2 -d1);
idxy = find(y<= L2/2 +d2 & y>=L2/2 -d2);
idxz = find(z);
idx = intersect(idxx, idxy);
idx = intersect(idx,idxz);
[i, j, k] = ind2sub(size(x),idx);
hz(i,j,k) = -1;
end