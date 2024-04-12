function calculateInitialFRFs(degree, nsub, startFreq, stopFreq, ...
    discretization, sAlpha, sBeta, filename)
[mesh, material, optimization] = standardParameters;
mesh.degree = degree;
mesh.nsub = nsub;
material.alpha = sAlpha;
material.beta = sBeta;
problem_data = cell(11,1);
for i=1:11
    problem_data{i} = square_shell_problem(i, [1 2 3 4]);
end

% Matrices
[K,M,C,dr_dofs, free_dofs, dr_values] = parseInitialMatrices(problem_data{1},mesh,material, optimization);
% Vectors
F = cell(11,1);
for i=1:11
    F{i} = parseInitialForce(problem_data{i},mesh,material);
end

% Boundary data


% Spectral calculations
frequency_array = linspace(startFreq,stopFreq,(stopFreq-startFreq)/discretization +1);

AIP = cell(11,1);
v2 = cell(11,1);
tmp = zeros(length(frequency_array),1);
tmp2 = tmp;
for i=1:11
    for j=1:length(frequency_array)
        f = frequency_array(j);
        omega = 2*pi*f;
        Kd = K +1j*omega*C -omega*omega*M;
        u = SolveDirichletSystem(Kd,F{i}, dr_dofs, free_dofs, dr_values);
        tmp(j) = real(0.5*omega*omega*u'*C*u);
        velocity = 1j*omega*u;
        tmp2(j) = real(velocity'*velocity);
    end
    AIP{i} = tmp;
    v2{i} = tmp2;
end
save(filename+string('.mat'),'frequency_array','AIP','v2','F');      
end