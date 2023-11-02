function dm = PlateMassSensitivities(lambda, u, t, msh, sp, lm, Me)
m_shape = [sp.nsh_max, sp.nsh_max];
dk = zeros(msh.nel,1);
for i=1:msh.nel
    dofs = lm(i,:);
    ell = lambda(dofs).';
    you = u(dofs);
    dm = reshape(Me(i,:),m_shape);
    dk(i) = ell*(dm)*you;
end

end