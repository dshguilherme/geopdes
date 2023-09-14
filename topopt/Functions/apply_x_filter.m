function xPhys = apply_x_filter(filter_options, x)
h = filter_options.h;
Hs = filter_options.Hs;
shape = filter_options.subshape;
switch filter_options.type
    case "simple"
      xPhys = x;
        
    case "density"
        xPhys = conv2(reshape(x,shape),h,'same')./Hs;
        xPhys = xPhys(:);
end
end