function [K,M,C, dr_dofs, free_dofs, dr_values] = parseInitialMatrices(problem, mesh, material)
degree = mesh.degree;
nsub = mesh.nsub;
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build spaces
[geometry, msh, sp] = buildSpaces(problem,method_data);

% Matrices
msh1 = msh_precompute(msh);
sp1 = sp_precompute_param(sp,msh1, 'value', true,'gradient', true, 'hessian', true);
for idim = 1:msh1.rdim
      x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
end
% Stiffness Elementary Matrix
[Bke, Ske] = op_KL_shells_elements(sp1,sp1,msh1,problem.E_coeff(x{:}), problem.nu_coeff(x{:}));
% Mass Elementary Matrix
Me = op_u_v_elements(sp1, sp1, msh1, material.rho);
% Location matrix (connectivity)
lm = sp1.connectivity';

% Stiffness, mass and damping
t = optimization.thickness;
[K, M] = shellMatricesFromElements(Bke, Ske, Me, lm, t, material.young, material.rho);
C = material.alpha*M +material.beta*K;

dr_dofs = [];
for iside = 1:numel(problem.drchlt_sides)
  side = problem.drchlt_sides(iside);
  components = problem.drchlt_components{iside};
  for icomp = components
    dr_dofs = union (dr_dofs, sp.boundary(side).dofs(sp.boundary(side).comp_dofs{icomp}));
  end
end

free_dofs = setdiff(1:sp.ndof, dr_dofs);
end