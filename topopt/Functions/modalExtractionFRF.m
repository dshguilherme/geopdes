function FRF = modalExtractionFRF(frequency_array, amount, t_val)
tval = t_val;
load('init_shell.mat')
tval = (tmax-tmin)*tval/100 +tmin;
t_val = apply_x_filter(filter_options,tval);
Ks = shellStiffnessFromElements(Bke, Ske, lm, t_val, t_val, YOUNG, modo);
M = shellMassFromElements(Me, lm, t_val, t_val, RHO, modo);
alpha_ = 1e-5;
beta_ = 1e-4;

FRF = zeros(numel(frequency_array),4);
receptance = extractReceptance(Ks,M,dr_dofs,amount, alpha_, beta_);
for i=1:length(frequency_array)
    f = frequency_array(i);
    omega = 2*pi*f;
    Aij = receptance(omega);
    X = Aij*F;
    V = 1j*omega*X;
    FRF(i,1) = V'*V;
    FRF(i,2) = 0.5*F'*V;
    FRF(i,3) = abs(F'*X);
    FRF(i,4) = real(V'*blkdiag(R0,R0,R0)*V);
end
end