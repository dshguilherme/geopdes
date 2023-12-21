function [f0val, fval] = EvalMultiShellObjectives(xval)
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
switch objective_function
    case "v2_rms"
        f0val = mean(V2_scaled);
end
end