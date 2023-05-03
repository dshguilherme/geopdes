function f = forceHeat(x,y,ind)
d = 0.1;
switch ind
        case 1
            f = zeros(1, size(x,1), size(x,2));
            idxx = 1:size(x,1);
            idxy = find(y<= 0.5 && y>= 0.5-d);
            idx = intersect(idxx,idxy);
            [i, j] = ind2sub(size(x),idx);
            f(1,i,j) = 1;
    end
end