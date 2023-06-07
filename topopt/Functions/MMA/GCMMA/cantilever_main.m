clearvars
clc

degree = 1;
nsub = [240 80];
vol_frac = 0.49;
freq = 100;
omega = 2*pi*freq;
rmin = 3;
change_min = 1e-4;
iter_max = 300;

YOUNG = 210e9;
RHO = 7860;
alpha = 0;
beta = 0.1/omega;
eta = 0.9;

[geometry, msh, sp, Ke, Me, Ve, F, lm] = initialize_cantilever(degree, nsub);
[h, Hs] = density_filter(rmin, nsub);

gccantileverinit;

% GCMMA loop
change = 0.5; loop = 0;
while change > change_min && loop < iter_max
    loop = loop+1;
    gccantilever;
    change = max(abs(xval-xold1));
    x_plot = conv2(reshape(xmma,nsub),h,'same')./Hs;
    volum = mean(xmma);
    fprintf(' Iteration.:%5i | Objective.:%11.2f | Vol.:%7.3f | Change.:%7.3f\n', ...
    loop, f0val, volum,change);
    colormap(gray); imagesc(1-rot90(x_plot)); caxis([0 1]); axis equal; axis off; drawnow;
end



