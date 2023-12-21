function problem_data = hemispherical_shell_problem(R,theta)
clear problem_data
crv = nrbcirc(R,[0 0 0],0,pi/2 -deg2rad(theta));
srf = nrbrevolve(crv,[0 0 0],[0 -1 0],pi/2);
srf = nrbtform(srf, vecrotz(pi/2));
srf = nrbtform(srf, vecrotx(pi/2));
srf = nrbtform(srf, vecroty(pi/2));

problem_data.geo_name = srf;
problem_data.drchlt_sides = [3];
problem_data.drchlt_components = {[3]};

% Physical parameters
E = 1;
nu = 0.3;
thickness = 1;
density = 1;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;
problem_data.density = density;

% Source and boundary terms
hx = @(x, y, z) zeros(size(x));
hy = @(x, y, z) zeros(size(x));
hz = @(x, y, z) -1*ones(size(x));

problem_data.f       = @(x, y, z, ind) cat(1, ...
    reshape (hx (x,y,z), [1, size(x)]), ...
    reshape (hy (x,y,z), [1, size(x)]), ...
    reshape (hz (x,y,z), [1, size(x)]));

end