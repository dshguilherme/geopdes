function problem_data = compliant_mechanism(L, h)
L = 1;
hh = L/2;
d = L/100;
k_in = 1;
k_out = k_in/100;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);
% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
problem_data.symm_sides   = [4];
% Physical parameters
E  =  1; nu = .3; Emin = 1e-9;
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));
problem_data.rho = @(x, y) ones(size(x));
% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));


end