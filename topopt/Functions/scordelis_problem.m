function problem_data = scordelis_problem
clear problem_data
radius = 25.0;
theta = 2.0 / 9.0 * pi;

nrb = nrbreverse (nrbcirc(radius, [0 0 0], pi/2 - theta, pi/2 + theta));
srf = nrbextrude (nrbtform(nrb,vecrotx(pi/2)), [0 50 0]);
problem_data.geo_name = srf;

% Type of boundary conditions for each side of the domain
% Only homogeneous Dirichlet conditions have been implemented so far.
problem_data.drchlt_sides = [3 4];
problem_data.drchlt_components = {[1 3] [1 3]};

% Physical parameters
E = 1;
nu = 0.0;
thickness = 1;
density = 1;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;
problem_data.density = density;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -90*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

end