function problem_data = square_shell_problem
L1 = 1; %pi;
L2 = 1; %exp(1);
p11 = [0 0 1]; p12 = [L1 0 1]; p21 = [0 L2 1]; p22 = [L1 L2 1];
problem_data.geo_name = nrb4surf (p11,p12,p21,p22);

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [1 2 3 4];
problem_data.drchlt_components = {[1 2 3] [1 2 3] [1 2 3] [1 2 3]};

% Physical parameters
E = 1;
nu = 0.3;
rho = 1;
thickness = 1;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;
problem_data.density = rho;
problem_data.rho_coeff = @(x,y,z) rho*ones(size(x));
% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @forceSquarePlateCentered;

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));
end