function xval = OC(f1, init_mat, filter_options)
    load(init_mat);
    loop = 0;
    change = 1;    
 while loop <= iter_max && change > change_min
        loop = loop +1;
        [f0val, df0dx, fval, dfdx] = f1(xval);
        ell1 = 0; ell2 = 1e9; move = 0.2;
        x = xval;
        dc = df0dx;
        dv = dfdx;
    while(ell2-ell1)/(ell1+ell2) > 1e-3
        mid = 0.5*(ell2+ell1);
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc(:)./dv(:)/mid)))));
        % Filter new value
        xPhys = apply_x_filter(filter_options,xnew);
        % Check restriction
        if sum(xPhys(:)) > vol_frac*length(x)
            ell1 = mid;
        else
            ell2 = mid;
        end
    end
    change = max(abs(xnew(:)-x(:)));
    xval = xnew;
    x_plot = reshape(xPhys,nsub);
    fprintf(' Iteration: %3i | Objective:%10.4f | Volume: %4.2f | Change:%7.4f\n', ...
    loop, f0val, mean(xPhys), change);
    colormap(gray); imagesc(1-rot90(x_plot)); caxis([0 1]); axis equal;axis off;drawnow;
 end
end