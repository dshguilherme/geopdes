function dm = PlateMassSensiSIMP (lambda, u, t, msh, sp, lm, Me, RHO, RHO_MIN)
m_shape = [sp.nsh_max, sp.nsh_max];
dm = zeros(msh.nel,1);
for i=1:msh.nel
    dofs = lm(i,:)';
    ell = lambda(dofs).';
    you = u(dofs);
    m_e = reshape(Me(i,:),m_shape);
    if t(i) > 0.1
        dm(i) = (RHO-RHO_MIN)*ell*m_e*you;
    elseif t(i) <= 0.1
        dm(i) =  9*(t(i)^8)*(RHO-RHO_MIN)*ell*m_e*you;
    end
end
end