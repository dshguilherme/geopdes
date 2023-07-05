function de = MassSensitivities(lambda, u, x, nel, lm, Me, RHO, RHO_MIN)
de = zeros(nel,1);
for e=1:nel
    dofs = lm(e,:)';
    m_e = squeeze(Me(e,:,:));
    if x(e) > 0.1
        de(e) = (RHO-RHO_MIN)*(lambda(dofs)')*m_e*u(dofs);
    elseif x(e) <= 0.1
        de(e) = 9*(x(e)^8)*(RHO-RHO_MIN)*(lambda(dofs)')*m_e*u(dofs);
    end
end
end