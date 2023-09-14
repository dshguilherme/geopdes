function initialize_biharmonic_problem(parameters, problem_data);
data_names = fieldnames(parameters);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= parameters.(data_names{iopt});']);
end
clear method_data
method_data.degree     = [degree degree];     % Degree of the bsplines
method_data.regularity = [degree-1 degree-1];     % Regularity of the splines
method_data.nsub       = nsub;     % Number of subdivisions
method_data.nquad      = [degree+1 degree+1];     % Points for the Gaussian quadrature rule

% Build Spaces
[geometry, msh, sp] = buildKirchoffSpace(problem_data,method_data);

% Pre-calculate element volume, stiffness, mass and dofs
msh1 = msh_precompute(msh);
sp1 = sp_precompute(sp,msh1, 'value', true,'gradient', true, 'hessian', true);
   for idim = 1:msh1.rdim
      x{idim} = reshape (msh1.geo_map(idim,:,:), msh1.nqn, msh1.nel);
    end
Ke = op_gradgradu_gradgradv_elements(sp1,sp1,msh1,problem_data.c_diff(x{:}));
Me = op_u_v_elements(sp1, sp1, msh1, problem_data.c_diff(x{:}));
lm = sp1.connectivity';

end