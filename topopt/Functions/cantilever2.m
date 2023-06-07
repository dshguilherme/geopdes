function [f0val, df0dx, fval, dfdx] = cantilever2(x, msh, sp, Ke, Me, F, ...
    Ve, lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0, h, Hs, eta)
[f0val, fval, u, uu] = cantilever1(x, msh, sp, Ke, Me, F, Ve, ...
    lm, YOUNG, RHO, omega, alpha, beta, vol_frac, W0, eta);
df0dx = zeros(msh.nel,1);
df1dx = df0dx;
for e=1:msh.nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    m_e = squeeze(Me(e,:,:));
    
    dk = (3*x(e)^2)*(YOUNG -1e-3)*k_e;
    if x(e) > 0.1
        dm = (RHO -1e-3)*m_e;
    elseif x(e) <= 0.1
        dm = 9*(x(e)^8)*(RHO -1e-3)*m_e;
    end
    dc = alpha*dm +beta*dk;
    
    dkd = (dk +1j*dc -omega*omega*dm);
    
    df0dx(e) = -0.5*omega*real(1i*(u(dofs)')*dkd*u(dofs));
    df1dx(e) = -(uu(dofs)'*dk*uu(dofs));
end
% Filtering chain-rule
nsub = msh.nel_dir;
df0dx = 10/(log(10)*f0val)*df0dx;
df0dx = reshape(eta.*df0dx+(1-eta).*df1dx, nsub);
df0dx = conv2(df0dx./Hs,h,'same');
df0dx = df0dx(:);

dfdx = 100*Ve'/sum(Ve);
dfdx = reshape(dfdx,nsub);
dfdx = conv2(dfdx./Hs,h,'same');
dfdx = dfdx(:)';
end