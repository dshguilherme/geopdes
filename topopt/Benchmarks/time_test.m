Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
Me = Ke;
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, m, dofs] = elementaryMatrices(i, sp, sp, msh, ...
        problem_data.lambda_lame, problem_data.mu_lame, problem_data.rho);
    Ke(i,:,:) = k;
    Me(i,:,:) = m;
    lm(i,:) = dofs;
end
[Ks, C, M] = SIMPMatrices(sp, msh, lm, Ke, Me, alpha_, beta_,YOUNG, ...
    YOUNG_MIN, RHO, RHO_MIN, x);

lambda_lame = problem_data.lambda_lame;
mu_lame = problem_data.mu_lame;

mat = op_su_ev_tp (sp, sp, msh, lambda_lame, mu_lame);