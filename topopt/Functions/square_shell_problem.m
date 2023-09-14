function problem_data = square_shell_problem
base = 10;
p11 = [0 0 1]; p12 = [base 0 1]; p21 = [0 base 1]; p22 = [base base 1];
problem_data.geo_name = nrb4surf (p11,p12,p21,p22);

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [1 2 3 4];
problem_data.drchlt_components = {[2 3] [2 3] [1 3] [1 3]};

% Physical parameters
E = 1e6;
nu = 0.;
thickness = (12 * (1 - nu*nu) / E)^(1/3);

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) 4*16*pi^4/(base^4)*sin(2*pi*x/base).*sin(2*pi*y/base);

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));
end