clearvars
clc
close all

%% Spline Space
L = 1;
h = 0.5;
vol_frac = 0.49; % Initial mass distribution
degree = 1; % Degree of the base functions you want
nsub = [150 50]; % Number of sub-divisions

regularity = [degree-1 degree-1];
geo_name = nrb4surf([0 0 vol_frac], [L 0 vol_frac], [0 h vol_frac], [L h vol_frac]);
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
nurbs = nrbkntins (nurbs, nknots);
density_space = geo_load (nurbs);

%% IGA Geometry
problem_data = cantilever_beam_special(L,h,nurbs);
clearvars -except problem_data
% Mesh Parameters
degree = 1;
parameters.degree = degree;
nsub = [90 30];
nsub = [150 50];


clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

% Pre-calculate element volume, stiffness, mass and dofs
K = op_su_ev_tp(sp, sp, msh, problem_data.lambda_lame, problem_data.mu_lame);

F = buildForce(sp, msh, problem_data);
F = 100000*F/sum(F);

[free_dofs, dr_dofs] = grab_cantilever_dofs(sp);

dr_values = zeros(length(dr_dofs),1);
us = SolveDirichletSystem(K,F,dr_dofs,free_dofs,dr_values);