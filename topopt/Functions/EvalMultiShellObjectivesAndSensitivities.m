function [f0val, df0dx, fval, dfdx] = EvalMultiShellObjectivesAndSensitivities(xval)
%% Filtering
tval =  xval;
load('init_shell_multi.mat');
tval = (tmax-tmin)*tval/100 +tmin;
t = apply_x_filter(filter_options,tval);
%% Assembly
Ks = shellStiffnessFromElements(Bke, Ske, lm, t, t, YOUNG, modo);
M = shellMassFromElements(Me, lm, t, t, RHO, modo);
C =  @(omega) alpha_(omega)*M +beta_(omega)*Ks;
Kd = @(omega) Ks +1j*omega*C(omega) -omega*omega*M;

%% Solve problems
dr_values = zeros(length(dr_dofs),1); 
u = cell(size(freq));
velocity = u;
for i=1:length(u)
    u{i} = SolveDirichletSystem(Kd(omega(i)), F, dr_dofs, free_dofs, dr_values);
    velocity{i} = -1j*omega(i)*u{i};
end

%% Objective Functions and Constraints
V2 = zeros(size(freq));
for i=1:3
    bigR = blkdiag(R{i},R{i},R{i});
    V2(i) = real(velocity{i}'*bigR*velocity{i});
end
V2_scaled = 100*V2./V0';

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
dc = @(omega) alpha_(omega)*dm +beta_(omega)*dk;

% Restrictions
dh1dt = 100*(RHO.*Ve')/M_max;
dh2dt = 100*(-RHO.*Ve')/M_min;

switch objective_function
    case "v2_rms"
        f0val = mean(V2_scaled);
        df0dx = zeros(prod(nsub),1);
        % Adjoint Problems
        for i=1:3
            you = u{i};
            lhs = -2*omega(i)*omega(i)*(you')*(bigR);
            ell = Kd(omega(i))\lhs.';
            dkd = dk +1j*omega(i)*dc(omega(i)) -omega(i)*omega(i)*dm;
            df_you = CalculateSensivities(ell,you,lm,dkd);
            df0dx = df0dx +100*df_you/(3*V0(i));
        end
end
dfdx = [dh1dt; dh2dt]; 
[df0dx, dfdx] = apply_sensi_filter(filter_options, t, df0dx, dfdx);
df0dx = df0dx*dfactor;
end