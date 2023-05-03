function f =  forceForceInverter(x,y,ind)
d = 0.1;
f = zeros(1,size(x,1),size(x,2));
idxx = find (x<=d);
idxy = find(y<=0.5 & y>=0.4);
idx = intersect(idxx,idxy);
[i,j] = ind2sub(size(x),idx);
f(1,i,j) = 1;
end