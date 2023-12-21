function receptance = extractReceptance(Ks,M,dr_dofs,amount,alpha,beta)

free_dofs = setdiff(1:length(Ks),dr_dofs);
autovec = zeros(length(Ks),amount);
[V, W] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),amount,'smallestabs');

autovec(free_dofs,:) = V;


w_n = sqrt(diag(W));
csi = diag(alpha./(2*w_n) +beta.*w_n/2);
w_n = diag(w_n);

mid = @(omega) 1./(W -omega*omega +2j*csi.*omega.*w_n);

receptance = @(omega) autovec*mid(omega)*autovec.';
end

