function [f0val, fval] = vesselEvalShellObjectives(io, xval)
t = (io.tmax-io.tmin)*xval/100 +io.tmin;
[Ks, M] = shellMatricesFromElements(io.Bke, io.Ske, io.Me, io.lm, t, io.YOUNG, io.RHO);
Ks(io.strip_dofs, io.strip_dofs) = Ks(io.strip_dofs,io.strip_dofs) +io.bK(io.b_dofs,io.b_dofs);

C = io.alpha_*M +io.beta_*Ks;
Kd = Ks +1j*io.omega*C -io.omega*io.omega*M;

% Solve
u = zeros(length(Kd),1);
u(io.free_dofs) = Kd(io.free_dofs,io.free_dofs)\io.F(io.free_dofs);
u(io.ms(:,2)) = u(io.ms(:,1));

%% Objective Function
switch io.objective_function
    case "v2_rms"
        velocity = -1j*io.omega*u;
        V2_rms = real(velocity'*velocity);
        V0 = io.v2_init;
        V2_scaled = 100*V2_rms/V0;
        f0val = V2_scaled;
    case "AIP"
        aW0 = io.aW_init;
        aW = real(0.5*io.omega*io.omega*(u'*C*u));
        aW_scaled = 100*aW/aW0;
        f0val = aW_scaled;
end
% Restrictions
M_max = io.RHO*sum(io.Ve.*io.thickness)*(1+io.max_add);
M_min = io.RHO*sum(io.Ve.*io.thickness)*(1-io.max_take);
Mass = io.RHO*sum(io.Ve'.*t); % Current mass
h1 = 100*(Mass - M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;

fval = [h1; h2];
end