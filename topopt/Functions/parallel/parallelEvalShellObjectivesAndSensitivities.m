function [f0val, df0dx, fval, dfdx] = parallelEvalShellObjectivesAndSensitivities(io, xval)

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
C = io.alpha_*M +io.beta_*Ks;

% Dynamic stiffness
Kd = Ks +1j*io.omega*C -io.omega*io.omega*M;

%% Solve problems
dr_values = zeros(length(io.dr_dofs),1);
u = SolveDirichletSystem(Kd, io.F, io.dr_dofs, io.free_dofs, dr_values);

%% Objective Function

velocity = -1j*io.omega*u;
% V2_rms = real(velocity'*io.R0*velocity);
V2_rms = real(velocity'*velocity);
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

%% Sensitivities
nel_dof = size(io.lm,2);
dk = zeros(io.msh.nel, nel_dof, nel_dof);
dm = dk;
for i=1:io.msh.nel
    ks = reshape(io.Ske(i,:), [nel_dof, nel_dof]);
    kb = reshape(io.Bke(i,:), [nel_dof, nel_dof]);
    me = reshape(io.Me(i,:), [nel_dof, nel_dof]);
    dk(i,:,:) = io.YOUNG*(ks +3*t(i)*t(i)*kb);
    dm(i,:,:) = io.RHO*me;
end

dc = io.alpha_*dm +io.beta_*dk;
dkd = dk +1j*io.omega*dc -io.omega*io.omega*dm;

% Restrictions
dh1dt = 100*(io.RHO.*io.Ve')/M_max;
dh2dt = 100*(-io.RHO.*io.Ve')/M_min;

% Chain Rules
f0val = V2_db;
% Adjoint problem
% lhs = -2*io.omega*io.omega*(u')*io.R0;
lhs = -2*io.omega*io.omega*(u');
ell = Kd\lhs.';
df0dx = CalculateSensivities(ell,u,io.lm,dkd);
tmp = 100/V0_db;
tmp2 = 10/(log(10)*V2_rms);
df0dx = tmp*tmp2*df0dx;

% df0dx = 100*df0dx/io.V0;

dfdx = [dh1dt; dh2dt];
[df0dx, dfdx] = apply_sensi_filter(io.filter_options, t, df0dx, dfdx);
% Scale back to 0-100 variable
df0dx = df0dx*io.dfactor;
end