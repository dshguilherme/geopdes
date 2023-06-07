function [f0val, fval, u, uu] = cantilever1(x, msh, sp, Ke, Me, F, Ve, ...
    lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0, eta)
% Pre-alocating vectors/solutions
[free_dofs, dr_dofs] = grab_cantilever_dofs(sp);
u = zeros(sp.ndof,1);
K = zeros(sp.ndof);
M = K;
for e=1:msh.nel
    k_e = (1e-3 +(x(e)^3)*(YOUNG-1e-3))*squeeze(Ke(e,:,:));
    if x(e) > 0.1
        m_e = (1e-3 +(x(e))*(RHO -1e-3))*squeeze(Me(e,:,:));
    elseif x(e) <= 0.1
       m_e = (1e-3 +(x(e)^9)*(RHO -1e-3))*squeeze(Me(e,:,:));
    end
        idx = lm(e,:)';
        K(idx,idx) = K(idx,idx) +k_e;
        M(idx,idx) = M(idx,idx) +m_e;
end
K = sparse(K); M = sparse(M);
C = alpha*M +beta*K;
% Solving
Kd = K+1i*omega*C-(omega^2)*M;
u(dr_dofs) = 0;
uu = u;
u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
uu(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);
active_input_power = real(0.5*omega*omega*(u')*C*(u));
W_scaled = 100*(100 +10*log10(active_input_power))/W0; 
compl = F'*uu;
f0val = eta*W_scaled +(1-eta)*compl;
V_ = sum(x.*Ve);
fval = 100*(V_/sum(Ve) -vol_frac);
end