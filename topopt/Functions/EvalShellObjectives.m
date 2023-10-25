function [f0val, fval] = EvalShellObjectives(xval)

%% Filtering
t = xval;
load('init_shell.mat');
t = apply_x_filter(filter_options,t);

%% Assembly
Ks = shellStiffnessFromElements(Bke, Ske, lm, t);
M = shellMassFromElements(Me, lm, t);
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

% Mean quadratic velocity
V2_rms = 0;
Area = sum(Ve);
for i=1:msh.nel
    dofs = lm(i,:);
    V2_rms = V2_rms + Ve(i)*(u(dofs)'*u(dofs));
end
V2_rms = (omega*omega/Area)*V2_rms;

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
        f0val = [V2_db_scaled, Cs_scaled];
    case "Initial"
        f0val = [Cs, V2_rms];
    case "eigenmax"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),1,'sm');
        vals = diag(vals);
        f0val = -vals(1);
    case "eigengap"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        f0val = vals(1) - vals(2);      
   
end


end