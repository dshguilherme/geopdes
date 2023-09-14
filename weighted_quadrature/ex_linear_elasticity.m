clearvars
clc
close all

%% Initialize Parameters
% Geometry
L = 1;
h = 0.5;
problem_data = cantilever_beam(L,h);

% Material
YOUNG = 1; % Steel
POISSON = 0.3;
rho0 = 1;
coeff  = @(x, y, z) rho0*ones(size(x));

% Mesh
degree = 5;
nsub = [120 40];

method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

% Assembly
tempo = tic;
xPhys = ones(nsub);
K = Elasticity_WQ(msh, sp, geometry, YOUNG, POISSON, xPhys);
M1 = Mass_SIMP_WQ(msh, sp.scalar_spaces{1}, geometry, rho0, xPhys);
M = blkdiag(M1,M1);
fprintf('Tempo de Assembly Careca Solver: %.2e [s] \n', toc(tempo))
tempo = tic;
msh1 = msh_precompute(msh);
sp1 = sp_precompute(sp,msh, 'value', true,'gradient', true, 'divergence', true);
   for idim = 1:msh1.rdim
      x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
    end
ke = op_su_ev_elements(sp1,sp1,msh1,problem_data.lambda_lame(x{:}), problem_data.mu_lame(x{:}));
me = op_u_v_elements(sp1,sp1,msh1,coeff(x{:}));
fprintf('Tempo de Assembly Matrizes Elementares: %.2e [s] \n', toc(tempo))

tempo = tic;
[ke1, me1] = elementarySIMPMatrices(sp,msh);
fprintf('Tempo de Assembly Matrizes Elementares 2: %.2e [s] \n', toc(tempo));

tempo = tic;
K_geopdes = op_su_ev_tp(sp, sp, msh, problem_data.lambda_lame, problem_data.mu_lame);
M_geopdes = op_u_v_tp(sp, sp, msh, coeff);
fprintf('Tempo de Assembly GeoPDEs: %.2e [s] \n', toc(tempo))
fprintf('Stiffness Matrix difference norm: %.2e \n', norm(K-K_geopdes,'fro'));
fprintf('Mass Matrix difference norm: %.2e \n', norm(M-M_geopdes,'fro'));

% F = buildForce(sp, msh, problem_data);
% % F = sparse(-100000*F/sum(F));
% 
% % Solving
% [free_dofs, dr_dofs] = grab_cantilever_dofs(sp);
% 
% dr_values = zeros(length(dr_dofs),1);
% u = SolveDirichletSystem(K,F,dr_dofs,free_dofs,dr_values);
% u2 = SolveDirichletSystem(K_geopdes,F,dr_dofs,free_dofs,dr_values);
% F'*u
% F'*u2

