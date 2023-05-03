function [Compliance, dc, ve] = mbbCompliance(x, u,sp, msh, lambda, mu, V_frac)
    Compliance = zeros(msh.nel,1);
    dc = Compliance;
    ve = Compliance;
    for e=1:msh.nel
        [Ke, dofs] = elementaryStiffness(e, sp, sp, msh, lambda, mu);
        Compliance(e) = (x(e)^3)*u(dofs)'*Ke*u(dofs);
        dc(e) = 3*(x(e)^2)* (u(dofs)')*Ke*u(dofs);
    end
end