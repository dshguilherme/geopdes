% function problem_data = pinched_hemisphere(radius)
radius = 10;
problem_data.geo_name = sphere_octant(radius);

% Physical parameters
E = 6.825e7;
nu = 0.3;
thickness = 0.04;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness; % This will change in GCMMA

degree = 3;
nsub = [35 35];

method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% problem_data.drchlt_sides = [3];
% problem_data.drchlt_components = {[1 2 3] [1 2 3]};


[geometry, msh, sp] = buildSpaces(problem_data,method_data);

dofs_1 = sp.boundary(1).dofs;
dofs_2 = sp.boundary(2).dofs;
dofs_3 = sp.boundary(3).dofs;
dofs_4 = sp.boundary(4).dofs;

% dofs_13 = intersect(dofs_1,dofs_3);
% dofs_23 = intersect(dofs_2, dofs_3);
% dr_dofs = unique([dofs_13; dofs_23]);
dr_dofs = dofs_3;
free_dofs = setdiff(1:sp.ndof, dr_dofs);
dr_values = zeros(size(dr_dofs));

F_y = intersect(dofs_1, dofs_4);
F_y = F_y(2);
F_z = intersect(dofs_2, dofs_4);
F_z = F_z(3);


K = op_KL_shells_tp (sp, sp, msh, problem_data.E_coeff, problem_data.nu_coeff, thickness);
F = sparse(length(K),1);
F(F_y) = -1;
F(F_z) = 1;

u = SolveDirichletSystem(K,F, dr_dofs, free_dofs, dr_values);

output_file = strcat('pinched_deg',num2str(degree));
vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, sp, geometry, vtk_pts, output_file, 'Displacement')
% end