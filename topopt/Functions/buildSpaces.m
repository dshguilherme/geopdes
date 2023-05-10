function [geometry, msh, space] = buildSpaces(problem_data,method_data)
    % Extract the fields from the data structures into local variables
    data_names = fieldnames (problem_data);
    for iopt  = 1:numel (data_names)
      eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
    end
    data_names = fieldnames (method_data);
    for iopt  = 1:numel (data_names)
      eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
    end

    % Construct geometry structure
    geometry = geo_load (geo_name);
    degelev  = max (degree - (geometry.nurbs.order-1), 0);
    nurbs    = nrbdegelev (geometry.nurbs, degelev);
    [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

    nurbs = nrbkntins (nurbs, nknots);
    geometry = geo_load (nurbs);

    % Construct msh structure
    rule     = msh_gauss_nodes (nquad);
    [qn, qw] = msh_set_quad_nodes (zeta, rule);
    msh      = msh_cartesian (zeta, qn, qw, geometry);

    % Construct space structure
    space_scalar = sp_nurbs (nurbs, msh);
    scalar_spaces = repmat ({space_scalar}, 1, msh.rdim);
    space = sp_vector (scalar_spaces, msh);
    clear space_scalar scalar_spaces
end