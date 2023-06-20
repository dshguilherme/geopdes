function [Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha, beta, ... 
    YOUNG, YOUNG_MIN, RHO, RHO_MIN, x)
u = zeros(sp.ndof, 1);
Ks = zeros(sp.ndof);
M = Ks;
for e=1:msh.nel
        k_e = (YOUNG_MIN +(x(e)^3)*(YOUNG - YOUNG_MIN))*squeeze(Ke(e,:,:));
        if x(e) > 0.1
            m_e = (RHO_MIN +x(e)*(RHO -RHO_MIN))*squeeze(Me(e,:,:));
        elseif x(e) <= 0.1
            m_e = (RHO_MIN +(x(e)^9)*(RHO -RHO_MIN))*squeeze(Me(e,:,:));
        end
        idx = lm(e,:)';
        Ks(idx,idx) = Ks(idx,idx) +k_e;
        M(idx,idx) = M(idx,idx) +m_e;
end
    Ks = sparse(Ks); M = sparse(M);
    C = alpha*M +beta*Ks;   
end