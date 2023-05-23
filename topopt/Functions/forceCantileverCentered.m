function f = forceCantileverCentered(x,y)
L = 1;
hh = 0.5;
d = hh/10;
h1 = hh/2 +d/2;
h2 = hh/2 -d/2;

f = zeros(2, size(x,1), size(x,2));
idxx = find(x>= L-d);
idxy = find(y<=h1 & y>=h2);
idx = intersect(idxx, idxy);
[i, j] = ind2sub(size(x),idx);
f(2,i,j) = -1;
        
end