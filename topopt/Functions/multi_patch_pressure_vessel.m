function problem_data = multi_patch_pressure_vessel(A,B)
    clear problem_data
    e1 = nrbcirc(A, [0 0 0],0,pi/2);
    e1.coefs(1,:) = e1.coefs(1,:)*0.5;
    e1 = nrbtform(e1, vectrans([-0.2350+B/2, 0, 0]));
    rx = vecrot(pi, [0 1 0]);
    e2 = nrbtform(e1,rx);
    P1 = e1.coefs(:,3);
    P2 = e2.coefs(:,3);
    line = nrbline(P1,P2);

    problem_data = struct('geo_name', {},'E_coeff',{},'nu_coeff', ...
        {},'thickness',{},'density',{},'p',{});
    calota1 = nrbrevolve(e1, [0 0 0], [1 0 0], 2*pi);
    problem_data(1).geo_name = calota1;
    calota2 = nrbrevolve(e2, [0 0 0], [1 0 0], 2*pi);
    problem_data(2).geo_name = calota2;
    cilindro = nrbrevolve(line, [0 0 0], [1 0 0], 2*pi);
    problem_data(3).geo_name = cilindro;

    %Physical parameters
    E = 1;
    nu = 0.3;
    thickness = 1;
    density = 1;
    P = 1;
    
    for i=1:3
        problem_data(i).E_coeff = @(x, y, z) E * ones(size(x));
        problem_data(i).nu_coeff = @(x, y, z) nu * ones(size(x));
        problem_data(i).thickness = thickness;
        problem_data(i).density = density;
        problem_data(i).p = @(x,y,ind) P*ones(size(x));
    end

    
end
