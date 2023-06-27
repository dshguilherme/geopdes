function initialize_cantilever(parameters)
L = 1;
hh = 0.5;

data_names = fieldnames(parameters);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= parameters.(data_names{iopt});']);
end

problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [];

% Physical parameters
E  =  1; Emin = 1e-3; rho0=1;
problem_data.E = E;
problem_data.Emin = Emin;
nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @forceCantileverCentered;
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
% problem_data.rho = @(x, y) rho0*ones(2, size (x, 1), size (x, 2));
problem_data.rho = @(x, y) rho0*ones(size(x));
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildSpaces(problem_data,method_data);

% Pre-calculate element volume, stiffness, mass and dofs
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
% f1 = @mbb_functions_and_sensitivities;
% f2 = @mbb_functions;
f1 = @cantilever_functions_and_sensitivities;
f2 = @cantilever_functions;

tmp = objective_function;
save('init.mat');
objective_function = "Cs0";
save('init.mat');
[Cs0 ~] = f2(xval);
objective_function = "W0";
save('init.mat');
[W0 ~] = f2(xval);
objective_function = tmp;
save('init.mat');

end
