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
tic
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
toc
tmp_msh = msh_precompute(msh);
Ve = (tmp_msh.element_size.^2)';
F = buildForce(sp, msh, problem_data);
F = 100000*F/sum(F);

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

[free_dofs, dr_dofs] = grab_cantilever_dofs(sp);

[h, Hs] = density_filter(rmin, nsub);

filter_options.type = philter;
filter_options.h = h;
filter_options.Hs = Hs;
filter_options.shape = nsub;

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