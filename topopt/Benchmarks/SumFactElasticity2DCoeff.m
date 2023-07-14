function Cij = SumFactElasticity2DCoeff(xPhys,x,y)

rho = SplineInterp2D(nurbs,x,y);
E = 1e-3 +(rho.^3)*(210e9 -1e-3);
lambda = 0.3*E/((1+0.3)*(1-0.6));
mu = E/(2*(1.3));
Cij = zeros(2,2,numel(rho));
M = zeros(3,2,2);
M(:,:,1) = [1 0; 0 0; 0 1];
M(:,:,2) = [0 0; 0 1; 1 0];
for i=1:numel(rho)
    D = lambda(i)*[1 1 0; 1 1 0; 0 0 0] +mu(i)*[2 0 0; 0 2 0; 0 0 1];
    for j=1:2
        for k=1:2
            Cij(:,:,i) = Cij(:,:,i) + M(:,:,j)'*D*M(:,:,k);
        end
    end
end

end