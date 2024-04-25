function [f0val, df0dx, fval, dfdx] = vesselEvalShellObjectivesAndSensitivities(io, xval)
t = (io.tmax-io.tmin)*xval/100 +io.tmin;
[Ks, M] = shellMatricesFromElements(io.Bke, io.Ske, io.Me, io.lm, t, io.YOUNG, io.RHO);
for i=1:length(io.strip_dofs)
    sd = io.strip_dofs{i}(:);
    Ks(sd,sd) = Ks(sd,sd)+io.bK{i};
end

C = io.alpha_*M +io.beta_*Ks;
Kd = Ks +1j*io.omega*C -io.omega*io.omega*M;

% Sensitivities
nel_dof = size(io.lm,2);
dk = io.YOUNG*(io.Ske +3*t.*t.*io.Bke);
dk = reshape(dk, io.msh.nel, nel_dof, nel_dof);
dm = io.RHO.*io.Me;
dm = reshape(dm, io.msh.nel, nel_dof, nel_dof);

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
        
        % Chain Rules
        lhs = -2*io.omega*io.omega*u';
        ell = Kd\lhs.';
        dc = io.alpha_*dm +io.beta_*dk;
        dkd = dk+1j*io.omega*dc -io.omega*io.omega*dm;
        dfdx = CalculateSensivities(ell,u,io.lm,dkd);
        tmp = 100/V0;
        dfdx = tmp*dfdx;
        df0dx = dfdx;
    case "AIP"
        aW0 = io.aW_init;
        aW = real(0.5*io.omega*io.omega*(u'*C*u));
        aW_scaled = 100*aW/aW0;
        f0val = aW_scaled;
        
        % Chain Rules
        ell = -0.5*1j*io.omega*u;
        dc = io.alpha_*dm +io.beta_*dk;
        dkd = dk +1j*io.omega*dc -io.omega*io.omega*dm;
        dfdx = CalculateSensivities(ell,u,io.lm,dkd);
        tmp = 100/aW0;
        dfdx = tmp*dfdx;
        df0dx = dfdx;
end
% Restrictions
M_max = io.RHO*sum(io.Ve.*io.thickness)*(1+io.max_add);
M_min = io.RHO*sum(io.Ve.*io.thickness)*(1-io.max_take);
Mass = io.RHO*sum(io.Ve'.*t); % Current mass
h1 = 100*(Mass - M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;
dh1dt = 100*(io.RHO.*io.Ve')/M_max;
dh2dt = 100*(-io.RHO.*io.Ve')/M_min;


fval = [h1; h2];
dfdx = [dh1dt'; dh2dt'];
df0dx = df0dx*io.dfactor;
end