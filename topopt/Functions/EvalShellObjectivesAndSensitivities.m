function [f0val, df0dx, fval, dfdx] = EvalShellObjectivesAndSensitivities(xval)
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
% bigR = blkdiag(R0,R0,R0);
% V2_rms = real(velocity'*bigR*velocity);
V2_rms = real(velocity'*velocity);
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
aW_scaled = 100*aW/aW0;

% Mass Restriction
M_max = RHO*sum(Ve.*thickness)*(1+maximum_to_add); % Maximum added mass
M_min = RHO*sum(Ve.*thickness)*(1-maximum_to_take); % Minimum added mass
Mass = RHO*sum(Ve.*t); % Current mass
h1 = 100*(Mass -M_max)/M_max;
h2 = 100*(-Mass +M_min)/M_min;


fval = [h1; h2];

nel_dof = size(lm,2);
%% Sensitivities
dk = zeros(msh.nel,nel_dof,nel_dof);
dm = zeros(msh.nel,nel_dof,nel_dof);
switch modo
    case "Continuous"
        for i=1:msh.nel
            ks = reshape(Ske(i,:),[nel_dof, nel_dof]);
            kb = reshape(Bke(i,:),[nel_dof, nel_dof]);
            me = reshape(Me(i,:), [nel_dof, nel_dof]);
            dk(i,:,:) = YOUNG*(ks +3*t(i)*t(i)*kb);
            dm(i,:,:) = RHO*me;
        end

    case "SIMP"
        for i=1:msh.nel
            ks = thickness*reshape(Ske(i,:),[nel_dof, nel_dof]);
            kb = (thickness^3)*reshape(Bke(i,:),[nel_dof, nel_dof]);
            me = reshape(Me(i,:), [nel_dof, nel_dof]);
            tmp = ((1e-6 -1)*t(i) +1)^2;
            dk(i,:,:) = 1e-6*YOUNG*(ks +kb)*(1-1e-6)/tmp;
            dm(i,:,:) = RHO*me;
        end 
end
dc = alpha_*dm +beta_*dk;
dkd = dk +1j*omega*dc -omega*omega*dm;

% Restrictions
dh1dt = 100*(RHO.*Ve')/M_max;
dh2dt = 100*(-RHO.*Ve')/M_min;

% Chain Rules
% dc_scaled = 100*dcs/Cs0;





%% Output
% "compliance", "scaled compliance", "v2_rms", "v2_scaled", "v2_db",
% "mixed", "History" or "Initial" as options for objective_function

switch objective_function
    case "compliance"
        f0val = Cs_scaled;
        % Adjoint problem -> self-adjoint
        ell = -us;
        df0dx = CalculateSensivities(ell,us,lm,dk);
        df0dx = 100*df0dx/Cs0;
    case "dynamic compliance"
        f0val = Cd_scaled;
        % Adjoint problem -> self-adjoint
        ell = -u;
        df0dx = CalculateSensivities(ell,u,lm,dkd);
        df0dx = 100*df0dx/Cd0;
    case "AIP"
        f0val = aW_scaled;
        % Adjoint problem -> self-adjoint
%         ell = -0.5*1j*omega*u;
        lhs = -0.5*omega*omega*u'*(C+C.');
        ell = Kd\lhs.';
        df0dx = CalculateSensivities(ell,u,lm,dkd);
%         db_tax = 100*(10/(log(10)*aW))/aW0_db;
        df0dx = 100*df0dx/aW0;
    case "v2_rms"
        f0val = V2_scaled;
        % Adjoint problem
%         lhs = -2*omega*omega*(u')*eye(length(u));
        lhs = -2*omega*omega*(u')*(bigR);
        ell = Kd\lhs.';
        df0dx = CalculateSensivities(ell,u,lm,dkd);
        df0dx = 100*df0dx/V0;
    case "v2_mix"
        f0val = neta*V2_db +(1-neta)*Cs_scaled;
        lhs = -omega*omega*(u')*(R0+R0.');
        ell = Kd\lhs.';
        dv2 = CalculateSensivities(ell,u,lm,dkd);
        tmp2 = 10/(log(10)*V2_rms);
        ell = -us;
        dcs = CalculateSensivities(ell,us,lm,dk);
        df0dx = 100*(neta*dv2*tmp2/V0_db +(1-neta)*dcs/Cs0);
    case "v2_db"
        f0val = V2_db;
        lhs = -2*omega*omega*(u');
%         lhs = -omega*omega*(u')*(bigR+bigR.');
        ell = Kd\lhs.';
        df0dx = CalculateSensivities(ell,u,lm,dkd);
        tmp = 100/V0_db;
        tmp2 = 10/(log(10)*V2_rms);
        df0dx = tmp*tmp2*df0dx;
    case "AIP_mix"
        f0val = neta*aW_scaled +(1-neta)*Cs_scaled;
        lhs = -0.5*omega*omega*u'*(C+C.');
        ell = Kd\lhs.';
        daw = CalculateSensivities(ell,u,lm,dkd);
        db_tax = 100*(10/(log(10)*aW))/aW0;
        daw = daw*db_tax;
        ell = -us;
        dcs = CalculateSensivities(ell,u,lm,dkd);
        df0dx = neta*daw +100*((1-neta)*dcs/Cs0);
        
    case "eigengap"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        u1 = u;
        u2 = u;
        u1(free_dofs) = vec(:,1);
        u2(free_dofs) = vec(:,2);
        f0val = 10*(vals(1) - vals(2))/-eig0;
        de1 = dk -vals(1)*dm;
        de2 = dk -vals(2)*dm;
        dfdx1 = CalculateSensivities(u1,u1,lm,de1);
        dfdx2 = CalculateSensivities(u2,u2,lm,de2);
        df0dx = 10*(dfdx1-dfdx2)/-eig0;
    case "eigenmax"
        [vec, vals] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),2,'sm');
        vals = diag(vals);
        u1 = u;
        u1(free_dofs) = vec(:,1);
        f0val = -vals(1);
        de1 = -(dk -vals(1)*dm);
        df0dx = CalculateSensivities(u1,u1,lm,de1);
end
dfdx = [dh1dt; dh2dt]; 
[df0dx, dfdx] = apply_sensi_filter(filter_options, t, df0dx, dfdx);
df0dx = df0dx*dfactor;
end