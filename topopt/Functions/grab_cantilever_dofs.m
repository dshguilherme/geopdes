function [free_dofs, dr_dofs] = grab_cantilever_dofs(sp)
dr_dofs = sp.boundary(1).dofs;
free_dofs = setdiff(1:sp.ndof, dr_dofs);
end