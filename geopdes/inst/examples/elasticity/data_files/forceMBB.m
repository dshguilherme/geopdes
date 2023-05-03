function f =  forceMBB(x,y,ind)
L = 50;
d = 0.04*L;
h = 0.32*L;
    switch ind
        case 1
            f = zeros(2, size(x,1), size(x,2));
            idxx = find(x<=d);
            idxy = find(y>=h-d);
            idx = intersect(idxx, idxy);
            [i, j] = ind2sub(size(x),idx);
            f(2,i,j) = -1;
        case 4
            f = zeros(2, size(x,1), size(x,2));
            idxx = find(x<=d);
            idxy = find(y>=h-d);
            idx = intersect(idxx, idxy);
            [i, j] = ind2sub(size(x),idx);
            f(2,i,j) = -1;
    end
end