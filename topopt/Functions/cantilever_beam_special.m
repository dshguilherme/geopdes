function problem_data = cantilever_beam_special(L,h,nurbs)
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 h], [L h]);
% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [];

% Physical parameters
YOUNG = 210e9; Emin = 1e-3; rho0=1e-3; RHO = 7860;
rho = @(x,y) SplineSpace2DCoeff(nurbs,x,y);
density = @(x,y) rho0 + ((rho(x,y)>0.1).*rho(x,y) +(rho(x,y)<=0.1).*(rho(x,y).^9))*(RHO - rho0);
E = @(x,y) Emin +(rho(x,y).^3)*(YOUNG - Emin);
problem_data.E = E;
problem_data.Emin = Emin;
nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E(x,y))/((1+nu)*(1-2*nu)));
problem_data.mu_lame = @(x, y) (E(x,y)/(2*(1+nu)));

% Source and boundary terms
problem_data.f = @forceCantileverCentered;
problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
% problem_data.rho = @(x, y) rho0*ones(2, size (x, 1), size (x, 2));
problem_data.rho = density;
end