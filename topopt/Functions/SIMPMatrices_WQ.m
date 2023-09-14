function [Ks, C, M] = SIMPMatrices_WQ(msh, sp, geometry, YOUNG, POISSON, RHO, ...
    alpha, beta, xPhys)

    Ks = Elasticity_WQ(msh, sp, geometry, YOUNG, POISSON, xPhys);
    M1 = Mass_SIMP_WQ(msh, sp.scalar_spaces{1}, geometry, RHO, xPhys);
    M = blkdiag(M1,M1);
    C = alpha*M +beta*Ks;
end