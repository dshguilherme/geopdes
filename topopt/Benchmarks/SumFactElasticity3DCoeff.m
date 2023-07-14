function Cij = SumFactElasticity3DCoeff(xPhys,x,y,z)

rho = SplineInterp3D(xPhys,x,y,z);
E = 1e-3 +(rho.^3)*(210e9 -1e-3);
lambda = 0.3*E/((1+0.3)*(1-0.6));
mu = E/(2*(1.3));
Cij = zeros(3,3,numel(rho));
M = zeros(6,3,3);
M(:,:,1) = [1 0 0; 0 0 0; 0 0 0; 0 1 0; 0 0 1; 0 0 0];
M(:,:,2) = [0 0 0; 0 1 0; 0 0 0; 1 0 0; 0 0 0; 0 0 1];
M(:,:,3) =  [0 0 0; 0 0 0; 0 0 1; 0 0 0; 1 0 0; 0 1 0];
D = zeros(6,6);
for ii=1:numel(rho)
    D = [lambda(ii)*ones(3)+ 2*mu(ii)*eye(3), zeros(3);zeros(3), mu(ii)*eye(3)];
    for jj=1:3
        for kk=1:3
            Cij(:,:,ii) = Cij(:,:,ii) + M(:,:,jj)'*D*M(:,:,kk);
        end
    end
end

end