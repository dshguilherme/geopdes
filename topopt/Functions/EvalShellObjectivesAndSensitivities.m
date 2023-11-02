function [f0val, df0dx, fval, dfdx] = EvalShellObjectivesAndSensitivities(xval)
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

%% Sensitivities

% Restrictions
dh1dt = RHO.*Ve';
dh2dt = -RHO.*Ve';

% Chain Rules
% dc_scaled = 100*dcs/Cs0;





%% Output
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

switch objective_function
    case "compliance"
        f0val = Cs;
        df0dx = dcs;
    case "scaled compliance"
        f0val = Cs_scaled;
        df0dx = dc_scaled;
    case "v2_rms"
        f0val = V2_rms;
        ell = Kd\(2*omega*omega*(u'.')); % Adjoint system solve
        dk = PlateStiffnessSensiContinuous(ell, u, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(ell, u, t, msh, sp, lm, Me, RHO);
        dc = alpha_*dm +beta_*dk;
        df0dx = real(dk +1j*omega*dc -omega*omega*dm);
    case "v2_scaled"
        f0val = V2_scaled;
        ell = Kd\(2*omega*omega*(u'.')); % Adjoint system solve
        dk = PlateStiffnessSensiContinuous(ell, u, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(ell, u, t, msh, sp, lm, Me, RHO);
        dc = alpha_*dm +beta_*dk;
        df0dx = real(dk +1j*omega*dc -omega*omega*dm);
        df0dx = 100*df0dx/V0;
    case "v2_db"
        f0val = V2_db_scaled;
        ell = Kd\(2*omega*omega*(u'.')); % Adjoint system solve
        dk = PlateStiffnessSensiContinuous(ell, u, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(ell, u, t, msh, sp, lm, Me, RHO);
        dc = alpha_*dm +beta_*dk;
        df0dx = real(dk +1j*omega*dc -omega*omega*dm);
        dv_db_constant = 10/(log(10)*V2_rms);
        df0dx = 100*dv_db_constant*df0dx/V0_db;
    case "mixed"
        f0val = neta*V2_db_scaled +(1-neta)*Cs_scaled;
        df0dx = neta*dv_db +(1-neta)*dc_scaled;
    case "eigenmax"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),1,'sm');
        vals = diag(vals);
        new_vec = zeros(length(free_dofs)+length(dr_dofs),1);
        new_vec(free_dofs,1) = vec(:,1);
        vec1 = new_vec(:,1);
        f0val = -vals(1);
        dk = PlateStiffnessSensiContinuous(vec1, vec1, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(vec1, vec1, t, msh, sp, lm, Me, RHO);
        df0dx = -dk +vals(1)*dm;
    case "eigengap"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        new_vec = zeros(length(free_dofs)+length(dr_dofs),2);
        new_vec(free_dofs,1) = vec(:,1);
        new_vec(free_dofs,2) = vec(:,2);
        vec1 = new_vec(:,1);
        vec2 = new_vec(:,2);
        f0val = vals(1) - vals(2);
        dke1 = PlateStiffnessSensiContinuous(vec1, vec1, t, msh, sp, lm, Bke, Ske, YOUNG);
        dke2 = PlateStiffnessSensiContinuous(vec2, vec2, t, msh, sp, lm, Bke, Ske, YOUNG);
        dme1 = PlateMassSensiContinuous(vec1, vec1, t, msh, sp, lm, Me, RHO);
        dme2 = PlateMassSensiContinuous(vec2, vec2, t, msh, sp, lm, Me, RHO);
        df0dx = (dke1 -vals(1)*dme1) -(dke2 -vals(2)*dme2);
    case "AIP"
        f0val = aW;
        ell = -0.5*1j*omega*u;
        dk = PlateStiffnessSensiContinuous(ell, u, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(ell, u, t, msh, sp, lm, Me, RHO);
        dc = alpha_*dm +beta_*dk;
        df0dx = real(dk +1j*omega*dc -omega*omega*dm);
    case "dB AIP"
        f0val = aW_db_scaled;
        ell = -0.5*1j*omega*u;
        dk = PlateStiffnessSensiContinuous(ell, u, t, msh, sp, lm, Bke, Ske, YOUNG);
        dm = PlateMassSensiContinuous(ell, u, t, msh, sp, lm, Me, RHO);
        dc = alpha_*dm +beta_*dk;
        dW_db_constant = 10/(log(10)*aW);
        dkd = real(dk +1j*omega*dc -omega*omega*dm);
        dW_db = 100*dW_db_constant*dkd/aW0_db;
        df0dx = dW_db;
        
end
dfdx = [dh1dt; dh2dt]; 
end