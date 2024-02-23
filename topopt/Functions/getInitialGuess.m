function initialGuess = getInitialGuess(K, M, dr_dofs, n_eigs, frequency, ...
                                        space, geometry, nsub)
free_dofs = setdiff(1:length(K),dr_dofs);
[V, W] = eigs(K(free_dofs, free_dofs),M(free_dofs,free_dofs),n_eigs,'smallestabs');
W = diag(sqrt(W));

[~, closestIndex] = min(abs(frequency-2*pi*W));
eigenvector = V(:,closestIndex);
padded_vector = zeros(length(K),1);

padded_vector(free_dofs) = eigenvector;
[eu, F] = sp_eval(padded_vector, space, geometry, nsub);

[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
grafo = surf(X,Y,squeeze(eu(3,:,:)),'visible','off');
initialGuess = grafo.CData;
close
initialGuess = rescale(initialGuess,0,100);
end