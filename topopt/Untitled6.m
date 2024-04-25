for i=1:1000
freq = i;
omega = 2*pi*freq;
Kd = Ks+1j*omega*C -omega*omega*M;
u_init = zeros(length(Kd),1);
F(free_dofs) = F(free_dofs) -Kd(free_dofs, dr_dofs)*u_init(dr_dofs);
u_init(free_dofs) = Kd(free_dofs, free_dofs) \ F(free_dofs);
u_init(ms(:,2)) = u_init(ms(:,1)); % Master/slave displacements
aW_init = real(0.5*omega*omega*(u_init'*C*u_init));
velocity0 = -1j*omega*u_init;
v2_init = real(velocity0'*velocity0);
aW(i) = aW_init;
v2(i) = v2_init;
end

[V, W] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),100,'smallestabs');
u_init = zeros(length(Kd),1);
u_init(free_dofs) = V(:,1);
u_init(ms(:,2)) = u_init(ms(:,1));
u1 = u_init(1:4512);
u2 = u_init(4513:4513+4512);
u3 = u_init(end-4511:end);
d = geo_deform(u_init,sp,geometry);
d1 = geo_deform(u1,sp1,geo1);
d2 = geo_deform(u2,sp2,geo2);
d3 = geo_deform(u3,sp3,geo3);
nrbplot(d1.nurbs,[30 30])
hold on
nrbplot(d2.nurbs,[30 30])
nrbplot(d3.nurbs,[30 30])

[V, W] = eigs(Ks(free_dofs,free_dofs),M(free_dofs,free_dofs),100,'smallestabs');
u_init = zeros(length(Ks),1);
for i=1:10
    figure
u_init(free_dofs) = V(:,i);
u_init(ms(:,2)) = u_init(ms(:,1));
d = geo_deform(u_init,sp,geometry);
nrbkntplot(d.nurbs)
end
