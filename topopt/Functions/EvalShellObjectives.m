function [f0val, fval] = EvalShellObjectives(xval)

%% Filtering
t = xval;
load('init_shell.mat');
t = apply_x_filter(filter_options,t);

%% Assembly
Ks = shellStiffnessFromElements(Bke, Ske, lm, t, t, YOUNG, modo);
M = shellMassFromElements(Me, lm, t, t, RHO, modo);
C = alpha_*M +beta_*Ks;
Kd = Ks +1j*C -omega*omega*M;

%% Solve problems
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd, F, dr_dofs, free_dofs, dr_values);
us = SolveDirichletSystem(Ks, F, dr_dofs, free_dofs, dr_values);

%% Objective Functions and Constraints

% Compliance
Cs = F'*us; % Static Compliance
Cs_scaled = 100*Cs/Cs0; % Scaled Compliance

% Dynamic Compliance
Cd = F'*u;

% Active Input Power
aW = 0.5*omega*real(1j*F'*u); % Active Input Power
aW_db = 100 +10*log10(aW);
aW0_db = 100+10*log10(aW0);
aW_db_scaled = 100*aW_db/aW0_db;

% Mean quadratic velocity
velocity = 1j*omega*u;
V2_rms = velocity'*velocity;
V2_scaled = 100*V2_rms/V0; % Scaled quadratic velocity
V2_db = 100 +10*log10(V2_rms); % dB mean quadratic velocity
V0_db = 100 +10*log10(V0); % V0 dB
V2_db_scaled = 100*V2_db/V0_db; % Scaled dB mean quadratic velocity

M_max = RHO*sum(Ve.*thickness)*(1+maximum_to_add); % Maximum added mass
M_min = RHO*sum(Ve.*thickness)*(1-maximum_to_take); % Minimum added mass
Mass = RHO*sum(Ve.*t); % Current mass
h1 = Mass -M_max;
h2 = -Mass +M_min;

fval = [h1; h2];

%% Output
switch objective_function
    case "compliance"
        f0val = Cs;
    case "scaled compliance"
        f0val = Cs_scaled;
    case "v2_rms"
        f0val = V2_rms;
    case "v2_scaled"
        f0val = V2_scaled;
    case "v2_db"
        f0val = V2_db_scaled;
    case "mixed"
        f0val = neta*V2_db_scaled +(1-neta)*Cs_scaled;
    case "History"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        f0val = [V2_rms, aW, Cd, vals(1), vals(2)-vals(1)];
    case "Initial"
        f0val = [Cs, V2_rms, aW];
    case "eigenmax"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),1,'sm');
        vals = diag(vals);
        f0val = -vals(1);
    case "eigengap"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        f0val = vals(1) - vals(2);      
    case "AIP"
        f0val = aW;
    case "dB AIP"
        f0val = aW_db_scaled;
end


end