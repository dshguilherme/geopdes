function get_spectrum(f_range, x, lm, Ke, Me, F, alpha, beta, free_dofs, YOUNG, RHO)
K = zeros(max(max(lm)));
M = K;
u = zeros(size(F));
for e=1:length(x)
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
K = sparse(K); M = sparse(M); F = sparse(F);
C = alpha*M +beta*K;
dr_dofs = setdiff(1:length(x),free_dofs);
u(dr_dofs) = 0;


H = zeros(length(f_range),1);
for i=1:length(f_range)
    f = f_range(i);
    omega = 2*pi*f;
    Kd = K+ 1j*omega*C -omega*omega*M;
    u(free_dofs) = Kd(free_dofs,free_dofs)\F(free_dofs);
    H(i) = 0.5*omega*omega*(u')*C*u;
end

mag = abs(H);
figure(1)
semilogy(f_range,mag)
grid on
end