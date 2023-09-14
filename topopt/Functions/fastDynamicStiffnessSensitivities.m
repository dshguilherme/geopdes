function dkd = fastDynamicStiffnessSensitivities(ke,me,lm,lambda,u, omega, ...
    xPhys, YOUNG, RHO, alpha_, beta_)
[sz1, sz2] = size(ke);
aux_size = [sqrt(sz2), sqrt(sz2)];
x = xPhys(:);
SIMP_penalty_stiff = (3*x.^2)*(YOUNG -1e-3);
SIMP_penalty_mass = (x<=0.1).*9.*(x.^8)*(RHO-1e-3) +(x>0.1).*(RHO -1e-3);
ke = SIMP_penalty_stiff.*ke;
me = SIMP_penalty_mass.*me;
ce = alpha_*me +beta_*ke;
dke = ke +1j*omega*ce -omega*omega*me;
clearvars ke ce me
dkd = arrayfun(@(i) (lambda(lm(:,i)).')*reshape(dke(i,:),aux_size)*u(lm(:,i)),1:sz1);
dkd = real(dkd(:));
end