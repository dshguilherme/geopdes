function initialize_problem(parameters, problem_data)
data_names = fieldnames(parameters);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= parameters.(data_names{iopt});']);
end
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

% Pre-calculate element volume, stiffness, mass and dofs
% [ke, me, lm] = elementarySIMPMatrices(sp, msh);
% [sz1, sz2] = size(ke);
Ke = zeros(msh.nel,sp.nsh_max,sp.nsh_max);
Me = Ke;
lm = zeros(msh.nel,sp.nsh_max);
for i=1:msh.nel
    [k, m, dofs] = elementaryMatrices(i, sp, sp, msh, ...
        problem_data.lambda_lame, problem_data.mu_lame, problem_data.rho);
    Ke(i,:,:) = k;
    Me(i,:,:) = m;
    lm(i,:) = dofs;
end
tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';
F = buildForce(sp, msh, problem_data);
% F = 100000*F/sum(F);
F = zeros(sp.ndof,1);
LL = F;

m = 1; % Number of Restrictions
n = msh.nel; % Number of variables
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = vol_frac*eeen;
xold1   = xval;
xold2   = xval;
xmin    = 0*eeen;
xmax    = eeen;
low     = xmin;
upp     = xmax;
c       = 10000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
raa0    = 0.01;
raa     = 0.01*eeem;
raa0eps = 0.000001;
raaeps  = 0.000001*eeem;
outeriter = 0;
maxoutit  = iter_max;
kkttol  = 0;

%% Boundary dofs
% [free_dofs, dr_dofs] = grab_cantilever_dofs(sp);

    x_dofs = sp.comp_dofs{1}; % Index of x dofs
    y_dofs = sp.comp_dofs{2}; % Index of y dofs
    symm_dofs = findSymm(sp,msh,problem_data);
    % Boundary 1
    b1_dofs = sp.boundary(1).dofs;
    b1x = intersect(x_dofs,b1_dofs); % x dofs of boundary 1 
    b1y = intersect(y_dofs,b1_dofs); % y dofs of boundary 1

    dy_dof = 0.5/msh.nel_dir(2);
    ny1 = ceil(0.5/dy_dof); % # of dofs of the clamp

    drch_b1 = [b1x(1:ny1), b1y(1:ny1)]';

    % Boundary 4
    b4_dofs = sp.boundary(4).dofs;

    b4x = intersect(x_dofs, b4_dofs);

    k1_dof = intersect(b1x, b4x); % dof of spring 1 and f_in

    % Boundary 2
    b2_dofs = sp.boundary(2).dofs;
    b2x = intersect(x_dofs, b2_dofs);

    k2_dof = intersect(b2x, b4x); % dof of the spring 2 and u_out

    LL(k2_dof) = 1;
    F(k1_dof) = 1;
    
    dr_dofs = drch_b1;
    free_dofs = setdiff(1:sp.ndof, [dr_dofs; symm_dofs]);

[h, Hs] = density_filter(rmin, nsub);

filter_options.type = philter;
filter_options.h = h;
filter_options.Hs = Hs;
filter_options.subshape = nsub;

Cs0 = 1;
aW0 = 1;
tmp = objective_function;
objective_function = "Initial";

f1 = @EvalSIMPObjectivesAndSensitivities;
f2 = @EvalSIMPObjectives;
save('init.mat')
[f0, ~] = f2(xval);
Cs0 = f0(1);
aW0 = f0(2);
objective_function = tmp;
save('init.mat');
end