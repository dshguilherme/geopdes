function dk = PlateStiffnessSensiContinuous(lambda, u, t, msh, sp, lm, Bke, Ske, YOUNG)
m_shape = [sp.nsh_max, sp.nsh_max];
dk = zeros(msh.nel,1);
for i=1:msh.nel
    dofs = lm(i,:)';
    ell = lambda(dofs).';
    you = u(dofs);
    Bend = 3*t(i)*t(i)*reshape(Bke(i,:),m_shape);
    Stress = reshape(Ske(i,:),m_shape);
    dk(i) = ell*YOUNG*(Bend+Stress)*you;
end
end