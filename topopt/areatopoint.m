% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
L = 1;
h = 1;
d = L/10;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 h], [L h]);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Physical parameters
problem_data.c_diff  = @(x, y) ones(size(x));

% Source and boundary terms
problem_data.f = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) zeros(size(x));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data 
method_data.degree     = [3 3];       % Degree of the splines
method_data.regularity = [2 2];       % Regularity of the splines
method_data.nsub       = [8 8];       % Number of subdivisions
method_data.nquad      = [4 4];       % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, ~, K, F] = solve_laplace_iso_more_outputs (problem_data, method_data);

% Boundary Conditions
sp = space;

b_dofs = sp.boundary(1).dofs; % dofs of side 1 of the boundary
dx_dof = L./sp.ndof_dir;
n_drch_dof = ceil(d./dx_dof);

dx_dof = L./sp.ndof_dir;
n_drch_dof = ceil(d./dx_dof);
drchlt_dofs = [b_dofs(1:n_drch_dof(1,1))];
drchlt_dofs = unique(drchlt_dofs(:));


u_drchlt = zeros(numel(drchlt_dofs),1); % For this case, homogeneous BCs.
u = zeros(sp.ndof,1);
u(drchlt_dofs) = u_drchlt;
free_dofs = setdiff(1:sp.ndof, [drchlt_dofs; symm_dofs]);
F(free_dofs) = F(free_dofs) -K(free_dofs, drchlt_dofs)*u_drchlt;
u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'areatopoint_Deg3_Reg2_Sub8';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('The result is saved in the file %s \n \n', output_file);
sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

[eu, FF] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(FF(1,:,:)), squeeze(FF(2,:,:)));
surf(X,Y,eu)

Compliance = F'*u;
