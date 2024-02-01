function [connectivity, coordinates, element_vals, spatial_vals] = makeParaviewData(xval, nsub, grafo)
idx1 = nsub(1);
idx2 = nsub(2);
element_vals = xval;
spatial_vals = rot90(reshape(xval,nsub));
connectivity = zeros(prod(nsub),4);
for i=1:idx1
    connectivity(i,:) = [i, i+1, i+idx2+1, i+idx2+2];
end
for i=1:idx2-1
    idx = (1:idx1)';
    connectivity(idx +i*idx2,:) = connectivity(idx,:)+i*(idx2+1);
end

X = grafo.XData(:);
Y = grafo.YData(:);
Z = grafo.ZData(:);
coordinates = [X,Y,Z];
end