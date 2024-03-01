function data_struct = getSpectralData(io)

data_struct.Bke = io.Bke;
data_struct.Ske = io.Ske;
data_struct.Me = io.Me;
data_struct.F = io.F;
data_struct.freq = io.freq;
data_struct.dr_dofs = io.dr_dofs;
data_struct.free_dofs = io.free_dofs;
data_struct.dr_values = io.dr_values;
data_struct.tmax = io.tmax;
data_struct.tmin = io.tmin;
data_struct.lm = io.lm;
data_struct.YOUNG = io.YOUNG;
data_struct.RHO = io.RHO;

end