function [errl2 h_size] = plate_l2_error(u, uex, space, mesh, problem_data)
msh = msh_precompute(mesh);
sp = sp_precompute(space,msh,'value', true,'gradient', true, 'hessian', true, 'divergence', true);
valu = sp_eval_msh(u,sp,msh,{'value','stress'}, problem_data.lambda_lame, problem_data.mu_lame);

for idir = 1:msh.rdim
    x{idir} = reshape (msh.geo_map(idir,:,:), msh.nqn*msh.nel, 1);
end
    s_uexx = reshape((feval (uex, x{:})), msh.nqn,msh.nel);
    w = msh.quad_weights .* msh.jacdet;
    s = valu{2};
    s_xx = squeeze(s(1,1,:,:));
    errl2_elem = sum(sum((s_xx - s_uexx).^2,1).*w);
    errl2  = sqrt (sum (errl2_elem));
%     errl2_elem  = sqrt (errl2_elem);
    h_size = max(msh.element_size);
end