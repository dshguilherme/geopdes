
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
c       = 1000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
raa0    = 0.01;
raa     = 0.01*eeem;
raa0eps = 0.000001;
raaeps  = 0.000001*eeem;
outeriter = 0;
maxoutit  = 100;
kkttol  = 0;

[free_dofs, drchlt_dofs, symm_dofs] = grab_mbb_dofs(sp, msh, problem_data);
[h, Hs] = density_filter(rmin, nsub);
dr_dofs = unique(union(drchlt_dofs, symm_dofs));

