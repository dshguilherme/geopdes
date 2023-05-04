% EX_PLANE_STRAIN_SQUARE: solve the plane-strain problem on a square.

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
% Physical domain, defined as NURBS map given in a text file
L = 30;
hh = 10;
d = 0.04*L;

problem_data.geo_name = nrb4surf([0 0], [L 0], [0 hh], [L hh]);

% Type of boundary conditions
problem_data.nmnn_sides   = [1];
problem_data.press_sides  = [];
problem_data.drchlt_sides = [];
% problem_data.drchlt_components = {[ 2]};
problem_data.symm_sides   = [1];

% Physical parameters
E  =  1; Emin = 1e-3;
problem_data.E = E;
problem_data.Emin = Emin;
nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @forceMBB;
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 1];     % Degree of the bsplines
method_data.regularity = [2 0];     % Regularity of the splines
method_data.nsub       = [120 10];     % Number of subdivisions
method_data.nquad      = [4 2];     % Points for the Gaussian quadrature rule


% 3) Initialize TopOpt parameters
vol_frac = 0.5;
x = vol_frac*ones(prod(method_data.nsub),1);
density = x;
n = method_data.nsub;
xx = linspace(0,L,n(1));
yy = linspace(0,hh,n(2));


% Density filtering parameters
rmin = 3;
[dy, dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(method_data.nsub),h,'same');



%% TopOpt Initial Conditions
loop = 0;
change = 1;
xold1 = x;
xold2 = x;
low = ones(size(x));
upp = ones(size(x));
while change > 0.01
    % 3) CALL TO THE SOLVER
    loop = loop+1;
    [geometry, msh, space, ~, K, F, symm_dofs] = solve_linear_elasticity_more_outputs (problem_data, method_data, density);

    sp = space;

    % Boundary Conditions
    y_dofs = sp.comp_dofs{2}; % Index of y dofs
    by_dofs = sp.boundary(2).dofs; % dofs of side 2 of the boundary

    dx_dof = L./sp.ndof_dir;
    n_drch_dof = ceil(d./dx_dof);

    drchlt_dofs = intersect(y_dofs, by_dofs);
    drchlt_dofs = drchlt_dofs(n_drch_dof(2)); 

    u_drchlt = zeros(numel(drchlt_dofs),1); % For this case, homogeneous BCs.

    u = zeros(sp.ndof,1);
    u(drchlt_dofs) = u_drchlt;
    free_dofs = setdiff(1:sp.ndof, [drchlt_dofs; symm_dofs]);
    F(free_dofs) = F(free_dofs) -K(free_dofs, drchlt_dofs)*u_drchlt;
    u(free_dofs) = K(free_dofs, free_dofs)\F(free_dofs);

    compliance = F'*u;
    [c_e, dc, d2c, Ve, v_xe, Re] = mbbCompliance(density, u, sp, msh, problem_data.lambda_lame, ...
        problem_data.mu_lame, E, Emin);
    n = method_data.nsub;
    xPhys = reshape(density,n);
    % Density filter
    dc_xy = reshape(dc,n);
    dc = conv2(dc_xy./Hs,h,'same');
    dv = reshape(Ve,n);
    dv = conv2(dv./Hs,h,'same');

    % Calc Lagrange Multiplier and update x
    ell1 = 0; ell2 = 1e9; move = 0.2;
    while(ell2-ell1)/(ell1+ell2) > 1e-3
        mid = 0.5*(ell2+ell1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc(:)./dv(:)/mid)))));
        xPhys = conv2(reshape(xnew,n),h,'same')./Hs;
        if sum(xPhys(:)) > vol_frac*length(x)
            ell1 = mid;
        else
            ell2 = mid;
        end
    end
    change = max(abs(xnew(:)-x(:)));
    density = xnew(:);
    x = xnew(:);
    
    
    % Check if Volume restriction is respected
%     V = sum(Ve);
%     V_ = V*vol_frac;
%     restriction = sum(v_xe./V_) -1;
%     
%     


%     [x_new, low, upp] = updateMMA(x, xold1, xold2, loop, compliance, restriction, dc, d2c, dv, vol_frac, low, upp); 
%     density = x_new;
%     xold2 = xold1;
%     xold1 = x;
%     change = max(abs(xnew(:)-x(:)));
    fprintf(' Iteration.:%5i | Compliance.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
        loop, compliance, mean(xnew(:)),change);
colormap(gray); imagesc(xx,yy,1-rot90(xPhys)); caxis([0 1]); axis equal; axis off; drawnow;

end



% % 4) POST-PROCESSING. 
% % 4.1) Export to Paraview
% output_file = 'mbbbeam_L50L_d3_r2_sub9';
% 
% vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
% fprintf ('results being saved in: %s \n \n', output_file)
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
%     problem_data.lambda_lame, problem_data.mu_lame)
% 
% % 4.2) Plot in Matlab. Comparison with the exact solution.
% [eu, FF] = sp_eval (u, space, geometry, vtk_pts);
% [X, Y]  = deal (squeeze(FF(1,:,:)), squeeze(FF(2,:,:)));
% 
% % subplot (1,2,1)
% quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
% title ('Numerical solution'), axis equal tight

