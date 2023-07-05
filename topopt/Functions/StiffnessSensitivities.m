function de = StiffnessSensitivities(lambda, u, x, nel, lm, Ke, YOUNG, YOUNG_MIN)
de = zeros(nel,1);
for e=1:nel
    dofs = lm(e,:)';
    k_e = squeeze(Ke(e,:,:));
    de(e) = 3*(x(e)^2)*(YOUNG-YOUNG_MIN)*((lambda(dofs)')*k_e*u(dofs));
end
end