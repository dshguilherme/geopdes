function xval = MMA(f1,init_mat, filter_options)
    load(init_mat);
    loop = 0;
    change = 1;
    while loop <= iter_max && change > change_min
        loop = loop +1;
        [f0val, df0dx, fval, dfdx] = f1(xval);
        df0dx2 = zeros(size(df0dx));
        dfdx2 = zeros(size(dfdx));
        [xmma,~,~,~,~,~,~,~,~,low,upp] = mmasub(m, n, loop, xval, xmin, xmax, ...
            xold1, xold2, f0val, df0dx, df0dx2, fval, dfdx, dfdx2, low, upp, ...
            a0, a, c, d);
        
        % Filter new value
        xPhys = apply_x_filter(filter_options,xmma);
        
        change = max(abs(xmma(:)-xval(:)));
        xold2 = xold1;
        xold1 = xval;
        xval = xmma;
        % Output Stuff
        x_plot = reshape(xPhys,nsub);
        fprintf(' Iteration: %3i | Objective:%10.4f | Volume: %4.2f | Change:%7.4f\n', ...
      loop, f0val, mean(xPhys), change);
        colormap(gray); imagesc(1-rot90(x_plot)); caxis([0 1]); axis equal;axis off;drawnow;
    end
end