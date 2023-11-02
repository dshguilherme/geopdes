function u = SolveDirichletSystem(K, F, dr_dofs, free_dofs, dr_values)
u = zeros(length(K),1);
u(dr_dofs) = dr_values;
F(free_dofs) = F(free_dofs) -K(free_dofs, dr_dofs)*u(dr_dofs);
u(free_dofs) = K(free_dofs,free_dofs)\F(free_dofs);
end