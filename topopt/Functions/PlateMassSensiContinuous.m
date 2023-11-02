function dm = PlateMassSensiContinuous(lambda, u, t, msh, sp, lm, Me, RHO)
m_shape = [sp.nsh_max, sp.nsh_max];
dm = zeros(msh.nel,1);
for i=1:msh.nel
    dofs = lm(i,:);
    ell = lambda(dofs).';
    you = u(dofs);
    m = RHO*reshape(Me(i,:),m_shape);
    dm(i) = ell*(m)*you;
end

end