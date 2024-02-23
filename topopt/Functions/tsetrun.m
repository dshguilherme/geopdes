degree = 2;
nsub = [30 30];
freq = 2500;
Fmag = 1e6;
force_type = [2];
fixed_sides = [1 2 3 4];
initial_thickness = 0.001;
max_add = .5;
max_take = 0.5;
iter_max = 100;
obj_function = "AIP";

filename = 'test';

optimizeSquareShell(degree, nsub, freq, Fmag, force_type, obj_function, ...
    fixed_sides, initial_thickness, max_add, max_take, iter_max, filename)