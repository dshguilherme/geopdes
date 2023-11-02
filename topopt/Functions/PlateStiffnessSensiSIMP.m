function dk = PlateStiffnessSensiSIMP(lambda, u, t, msh, sp, lm, Bke, Ske, YOUNG)
m_shape = [sp.nsh_max, sp.nsh_max];
dk = zeros(msh.nel,1);
for i=1:msh.nel
    dofs = lm(i,:)';
    ell = lambda(dofs).';
    you = u(dofs);
    Bend = reshape(Bke(i,:),m_shape);
    Stress = reshape(Ske(i,:),m_shape);
    Evoid = YOUNG*1e-6;
    num = 1 - 1e-6;
    den = (-num*t(i) +1)^2;
    dYoung = (Evoid*num)/den;
    dk(i) = ell*dYoung*(Bend+Stress)*you;
end
end