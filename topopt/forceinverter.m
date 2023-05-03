% EX_PLANE_STRAIN_SQUARE: solve the plane-strain problem on a square.
% x_dofs = sp.comp_dofs{1};
% y_dofs = sp.comp_dofs{2};
% b_dofs = sp.boundary(1).dofs;
% bx_dofs = intersect(x_dofs, b_dofs);
% by_dofs = intersect(y_dofs, b_dofs);
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
L = 1;
h = 1;
d = L/100;
k_in = 1;
k_out = k_in/100;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 h], [L h]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
problem_data.symm_sides   = [4];

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the bsplines
method_data.regularity = [2 2];     % Regularity of the splines
method_data.nsub       = [9 9];     % Number of subdivisions
method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, ~, K, F, symm_dofs] = solve_linear_elasticity_more_outputs (problem_data, method_data);

% Boundary Conditions
x_dofs = sp.comp_dofs{1}; % Index of x dofs
y_dofs = sp.comp_dofs{2}; % Index of y dofs

b_dofs = sp.boundary(1).dofs; % dofs of side 1 of the boundary

bx_dofs = intersect(b_dofs, x_dofs);
by_dofs = intersect(b_dofs, y_dofs);
dx_dof = L./sp.ndof_dir;
n_drch_dof = ceil(d./dx_dof);
drchlt_dofs = [bx_dofs(1:n_drch_dof(1,1)); by_dofs(1:n_drch_dof(1,2))];
drchlt_dofs = unique(drchlt_dofs);

k1_dof = bx_dofs(end); % Where to apply k1 and F

b_dofs = sp.boundary(2).dofs; % dofs of side 2 of the boundary
bx_dofs = intersect(b_dofs, x_dofs);
k2_dof = bx_dofs(end); % Where to apply k2


K(k1_dof,k1_dof) = K(k1_dof,k1_dof) +k_in; % Adding spring to k1
K(k2_dof, k2_dof) = K(k2_dof,k2_dof) +k_out; % Adding spring to k2
LL = F;
LL(k2_dof) = 1;
F(k1_dof) = 1;


u_drchlt = zeros(numel(drchlt_dofs),1); % For this case, homogeneous BCs.

u = zeros(sp.ndof,1);
u(drchlt_dofs) = u_drchlt;
free_dofs = setdiff(1:sp.ndof, [drchlt_dofs; symm_dofs]);
F(free_dofs) = F(free_dofs) -K(free_dofs, drchlt_dofs)*u_drchlt;
u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);

% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'forceinverter_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% 4.2) Plot in Matlab. Comparison with the exact solution.
[eu, FF] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(FF(1,:,:)), squeeze(FF(2,:,:)));

% subplot (1,2,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution'), axis equal tight
% subplot (1,2,2)
% eu2 = problem_data.uex (X, Y);
% quiver (X, Y, squeeze(eu2(1,:,:)), squeeze(eu2(2,:,:)))
% title ('Exact solution'), axis equal tight


Compliance = LL'*u;

% error_l2 = sp_l2_error (space, msh, u, problem_data.uex)

%!demo
%! ex_plane_strain_square

%!test
%! problem_data.geo_name = nrb4surf([0 0], [1 0], [0 1], [1 1]);
%! problem_data.nmnn_sides   = [];
%! problem_data.press_sides  = [];
%! problem_data.drchlt_sides = [1 2 3 4];
%! problem_data.symm_sides   = [];
%! E  =  1; nu = .3; 
%! problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
%! problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));
%! fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
%!      (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
%! fy = fx;
%! problem_data.f = @(x, y) cat(1, ...
%!                 reshape (fx (x,y), [1, size(x)]), ...
%!                 reshape (fy (x,y), [1, size(x)]));
%! problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
%! uxex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
%! uyex = @(x,y) sin(2*pi*x).*(sin(2*pi*y));
%! problem_data.uex = @(x, y) cat(1, ...
%!                 reshape (uxex (x,y), [1, size(x)]), ...
%!                 reshape (uyex (x,y), [1, size(x)]));
%! method_data.degree     = [3 3];     % Degree of the bsplines
%! method_data.regularity = [2 2];     % Regularity of the splines
%! method_data.nsub       = [9 9];     % Number of subdivisions
%! method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
%! [geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);
%! error_l2 = sp_l2_error (space, msh, u, problem_data.uex);
%! assert (msh.nel, 81)
%! assert (space.ndof, 288)
%! assert (error_l2, 2.60376176743492e-04, 1e-16)
