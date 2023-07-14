L = 60;
hh = 20;
d = 0.05*L;
degree = 1;
nsub = [60 20];

problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [1];

% Physical parameters
YOUNG = 210e9; Emin = 1e-3; rho0=1e-3; RHO = 7860;
rho = @(x,y) SplineSpace2DCoeff(nurbs,x,y);
density = @(x,y) rho0 + ((rho(x,y)>0.1).*rho(x,y) +(rho(x,y)<=0.1).*(rho(x,y).^9))*(RHO - rho0);
E = @(x,y) Emin +(rho(x,y).^3)*(YOUNG - Emin);

% Source and boundary terms
problem_data.f = @forceMBBDistributed;
% problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
% problem_data.g = @forceMBB;
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.rho = density;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp, u] = solve_linear_elasticity(problem_data,method_data);