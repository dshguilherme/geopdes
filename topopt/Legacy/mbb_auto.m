clearvars
clc
nsub = {[60 20], [120 40], [150 50], [180 60]};
rmin = [2 3 4 4];
f = cell(4,1);
x = f;
iter = zeros(4,1);
for i=1:4
    ns = nsub{i};
    nelx = ns(1); nely = ns(2);
    [ff, xx, ell] = top71(nelx, nely, 0.5, 3, rmin(i), 1);
    f{i} = ff;
    x{i} = xx;
    iter(i) = ell;
end

degree = [1 2 3];
figa = cell(4,3);
xiga = figa;
iteriga = zeros(4,3);
for j=1:3
    for i=1:4
        [ff, xx, ell] = mbb_beam(degree(j), nsub{i}, 0.5, rmin(i), 'OC', 0.01, 400);
        figa{i,j} = ff;
        xiga{i,j} = xx;
        iteriga(i,j) = ell;
    end
end