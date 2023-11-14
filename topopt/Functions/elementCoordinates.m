function x = elementCoordinates(problem_data, method_data)
    % Extract the fields from the data structures into local variables
    data_names = fieldnames (problem_data);
    for iopt  = 1:numel (data_names)
      eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
    end
    data_names = fieldnames (method_data);
    for iopt  = 1:numel (data_names)
      eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
    end
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);
nquad = ones(1,length(geo_name.knots));
rule = msh_gauss_nodes(nquad);
[qn, qw] = msh_set_quad_nodes(zeta, rule);
msh = msh_cartesian(zeta, qn, qw, geometry);
msh = msh_precompute(msh);
x = cell(1,msh.rdim);
for idim = 1:msh.rdim
    x{idim} = reshape(msh.geo_map(idim,:,:), msh.nqn, msh.nel);
end
end