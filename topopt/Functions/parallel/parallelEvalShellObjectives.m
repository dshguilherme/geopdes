function [f0val, fval] = fastEvalShellObjectives(io, xval)

% Scale back thickness to real values
tval = (io.tmax-io.tmin)*xval/100 +io.tmin;
% Apply the appropriate filtering
t = apply_x_filter(io.filter_options,tval);

%% Assembly
% Build Stiffness and Mass from elementary values
% Ks = shellStiffnessFromElements(io.Bke, io.Ske, io.lm, t,t, io.YOUNG, io.modo);
% M = shellMassFromElements(io.Me, io.lm, t, t, io.RHO, io.modo);
[Ks, M] = shellMatricesFromElements(io.Bke, io.Ske, io.Me, io.lm, t, io.YOUNG, io.RHO);
% Rayleigh damping matrix


% Dynamic stiffness
u = cell(io.nfreq,1);
f0val = 0;
for i=1:io.nfreq
C = io.alpha_(i)*M +io.beta_(i)*Ks;
Kd = Ks +1j*io.omega(i)*C -io.omega(i)*io.omega(i)*M;

%% Solve problems
dr_values = zeros(length(io.dr_dofs),1);
u{i} = SolveDirichletSystem(Kd, io.F, io.dr_dofs, io.free_dofs, dr_values);

%% Objective Function

switch io.objective_function
    case "v2_rms"
    velocity = -1j*io.omega(i)*u{i};
    V2_rms = real(velocity'*velocity);
    % V2_rms = real(velocity'*io.R0*velocity);
    V0 = io.u_init{i,2};
    V2_scaled = 100*V2_rms/V0;
    % V_db = 100 +10*log10(V2_rms);
    % V0_db = 100+10*log10(V0);
    % V2_db = 100*V_db/V0_db;

    % Chain Rules
    f0val = f0val + V2_scaled/io.nfreq;
    case "AIP"
        aW0 = io.aW_init{i};
        aW = real(0.5*io.omega(i)*io.omega(i)*(u{i}'*C*u{i}));
        aW_scaled = 100*aW/aW0;
        f0val = f0val +aW_scaled/io.nfreq;
end
        
% f0val = f0val + V2_db/io.nfreq;
end

% Restrictions
M_max = io.RHO*sum(io.Ve.*io.thickness)*(1+io.maximum_to_add);
M_min = io.RHO*sum(io.Ve.*io.thickness)*(1-io.maximum_to_take);
Mass = io.RHO*sum(io.Ve.*t); % Current mass
h1 = 100*(Mass - M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;

fval = [h1; h2];
end