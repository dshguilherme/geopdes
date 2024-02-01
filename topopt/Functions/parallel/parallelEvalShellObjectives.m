function [f0val, fval] = fastEvalShellObjectives(io, xval)

% Scale back thickness to real values
tval = (io.tmax-io.tmin)*xval/100 +io.tmin;
% Apply the appropriate filtering
t = apply_x_filter(io.filter_options,tval);
t = gpuArray(t);

%% Assembly
% Build Stiffness and Mass from elementary values
Ks = shellStiffnessFromElements(io.Bke, io.Ske, io.lm, t,t, io.YOUNG, io.modo);
M = shellMassFromElements(io.Me, io.lm, t, t, io.RHO, io.modo);
% Rayleigh damping matrix
C = io.alpha_*M +io.beta_*Ks;

% Dynamic stiffness
Kd = Ks +1j*io.omega*C -io.omega*io.omega*M;

%% Solve problems
dr_values = zeros(length(io.dr_dofs),1);
u = SolveDirichletSystem(Kd, io.F, io.dr_dofs, io.free_dofs, dr_values);

%% Objective Function

velocity = -1j*io.omega*u;
V2_rms = real(velocity'*velocity);
% V2_rms = real(velocity'*io.R0*velocity);
V2_scaled = 100*V2_rms/io.V0;
V_db = 100 +10*log10(V2_rms);
V0_db = 100+10*log10(io.V0);
V2_db = 100*V_db/V0_db;


% Restrictions
M_max = io.RHO*sum(io.Ve.*io.thickness)*(1+io.maximum_to_add);
M_min = io.RHO*sum(io.Ve.*io.thickness)*(1-io.maximum_to_take);
Mass = io.RHO*sum(io.Ve.*t); % Current mass
h1 = 100*(Mass - M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;

fval = [h1; h2];
f0val = V2_db;
end