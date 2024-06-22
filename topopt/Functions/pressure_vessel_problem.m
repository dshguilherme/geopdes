function problem_data = pressure_vessel_problem(A,B)
clear problem_data
e1 = nrbcirc(A, [0 0 0],0,pi/2);
e1.coefs(1,:) = e1.coefs(1,:)*0.5;
e1 = nrbtform(e1, vectrans([-0.2350+B/2, 0, 0]));
rx = vecrot(pi, [0 1 0]);
e2 = nrbtform(e1,rx);
P1 = e1.coefs(:,3);
P2 = e2.coefs(:,3);
line = nrbline(P1,P2);
glue = nrbglue(e1,line);

domain = nrbglue(glue,e2);
domain.knots = domain.knots/max(domain.knots);
domain_interp = domain_interpolation(47,domain);

srf = nrbrevolve(domain_interp,[0 0 0], [1 0 0], 2*pi);
srf2 = nrbtform(srf, vecroty(pi));


problem_data.geo_name = srf;
problem_data.drchlt_sides = [1 2];
problem_data.drchlt_components = {[1 2 3], [1 2 3]};

%Physical parameters
E = 1;
nu = 0.3;
thickness = 1;
density = 1;

problem_data.E_coeff = @(x, y, z) E * ones(size(x));
problem_data.nu_coeff = @(x, y, z) nu * ones(size(x));
problem_data.thickness = thickness;
problem_data.density = density;

% Source terms
P = 1;
problem_data.p = @(x,y,ind) P*ones(size(x));

end