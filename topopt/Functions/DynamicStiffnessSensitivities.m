function dkd = DynamicStiffnessSensitivities(lambda, omega, u, x, nel, lm, ...
    Ke, Me, YOUNG, YOUNG_MIN, RHO, RHO_MIN, alpha_, beta_)
de = zeros(nel,1);
for e=1:nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    m_e = squeeze(Me(e,:,:));
    dk = 3*(x(e)^2)*(YOUNG-YOUNG_MIN)*k_e;
    if x(e) <= 0.1
        dm = 9*(x(e)^8)*(RHO-RHO_MIN)*m_e;
    elseif x(e) > 0.1
        dm = (RHO-RHO_MIN)*m_e;
    end
    dc = alpha_*dm +beta_*dk;
    dkd(e) = real((lambda(dofs)')*(dk +1j*omega*dc -omega*omega*dm)*u(dofs));
end
end