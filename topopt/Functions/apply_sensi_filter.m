function [df0dx, dfdx] = apply_sensi_filter(filter_options, xPhys, df0dx, dfdx)
h = filter_options.h;
Hs = filter_options.Hs;
shape = filter_options.shape;
switch filter_options.type
    case "simple"
        tmp = reshape(df0dx.*xPhys,shape);
        tmp2 = reshape(max(1e-3,xPhys),shape);
        df0dx = conv2(tmp,h,'same')./Hs./tmp2;
        df0dx = df0dx(:);
    case "density"
        df0dx = conv2(reshape(df0dx,shape)./Hs,h,'same');
        df0dx = df0dx(:);
        dfdx = conv2(reshape(dfdx,shape)./Hs,h,'same');
        dfdx = dfdx(:)';
end
end