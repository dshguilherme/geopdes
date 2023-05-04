function [Compliance, dc, d2c, Ve, v_xe, Re] = mbbCompliance(x, u,sp, msh, lambda, mu, E, Emin)
    Compliance = zeros(msh.nel,1);
    dc = Compliance;
    Ve = Compliance;
    for e=1:msh.nel
        [Ke, dofs] = elementaryStiffness(e, sp, sp, msh, lambda, mu);
        Compliance(e) = (Emin +x(e)^3*(E-Emin))*u(dofs)'*Ke*u(dofs);
        dc(e) = -3*(x(e)^2)*(E-Emin)* (u(dofs)')*Ke*u(dofs);
        d2c(e) = 6*(x(e))*(E-Emin)*(u(dofs)')*Ke*u(dofs);
        Ve(e) = elementVolume(e,sp,msh);
    end
    v_xe = x.*Ve;
    Re = dc./Ve;
end