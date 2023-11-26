function [f0val, fval] = EvalShellObjectives(xval)

%% Filtering
tval =  xval;
load('init_shell.mat');
tval = (tmax-tmin)*tval/100 +tmin;
t = apply_x_filter(filter_options,tval);

%% Assembly
Ks = shellStiffnessFromElements(Bke, Ske, lm, t, t, YOUNG, modo);
M = shellMassFromElements(Me, lm, t, t, RHO, modo);
C = alpha_*M +beta_*Ks;
Kd = Ks +1j*omega*C -omega*omega*M;

%% Solve problems
dr_values = zeros(length(dr_dofs),1);
u = SolveDirichletSystem(Kd, F, dr_dofs, free_dofs, dr_values);
us = SolveDirichletSystem(Ks, F, dr_dofs, free_dofs, dr_values);

%% Objective Functions and Constraints

% Static Compliance
Cs = F'*us;
Cs_scaled = 100*Cs/Cs0;

% Dynamic Compliance
Cd = abs(F'*u);
Cd_scaled = 100*Cd/Cd0;

% Quadratic Velocity
velocity = -1j*omega*u;
V2_rms = real(velocity'*blkdiag(R0,R0,R0)*velocity);
V2_scaled = 100*V2_rms/V0;
V_db = 100 +10*log10(V2_rms);
V0_db = 100 +10*log10(V0);
V2_db = 100*V_db/V0_db;


% Active Input Power
% aW = 0.5*omega*real(1j*F'*u);
% aW = 0.5*real(velocity'*C*velocity);
aW = real(0.5*omega*omega*(u'*C*u));
aW_db = 100 +10*log10(aW);
aW0_db = 100 +10*log10(aW0);
aW_scaled = 100*aW_db/aW0_db;
% aW_scaled = 100*aW/aW0;

% Mass Restriction
M_max = RHO*sum(Ve.*thickness)*(1+maximum_to_add); % Maximum added mass
M_min = RHO*sum(Ve.*thickness)*(1-maximum_to_take); % Minimum added mass
Mass = RHO*sum(Ve.*t); % Current mass
h1 = 100*(Mass -M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;


fval = [h1; h2];

%% Output
switch objective_function
    case "compliance"
        f0val = Cs_scaled;
    case "dynamic compliance"
        f0val = Cd_scaled;
    case "v2_rms"
        f0val = V2_scaled;
    case "v2_db"
        f0val = V2_db;
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
        f0val = 10*(vals(1) - vals(2))/-eig0;      
    case "AIP"
        f0val = aW_scaled;
    case "dB AIP"
        f0val = aW_db_scaled;
    case "v2_mix"
        f0val =  neta*V2_scaled +(1-neta)*Cs_scaled;
    case "AIP_mix"
        f0val = neta*aW_scaled +(1-neta)*Cs_scaled;
end


end