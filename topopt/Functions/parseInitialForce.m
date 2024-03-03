function [F] = parseInitialForce(problem, mesh, material)
degree = mesh.degree;
nsub = mesh.nsub;
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build spaces
[geometry, msh, sp] = buildSpaces(problem,method_data);

% Force vector
F = op_f_v_tp (sp, msh, problem.f);
F = sparse(F);

end