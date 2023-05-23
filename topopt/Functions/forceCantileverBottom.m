function f = forceCantileverBottom(x,y)
L = 1;
hh = 0.5;
d = hh/10;

f = zeros(2, size(x,1), size(x,2));
idxx = find(x>=L-d);
idxy = find(y<=d);
idx = intersect(idxx, idxy);
[i, j] = ind2sub(size(x),idx);
f(2,i,j) = -1;
        
end