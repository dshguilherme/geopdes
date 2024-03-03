function batchOptimizeKLShell(problem_data, parameters, filename)

io = firstStepParallelKLShell(parameters,problem_data);
[xval, fobj, fres, x_history] = fastGCMMA_noplots(io);

%% Save Paraview and Results
grafo = nrbplot(io.geometry.nurbs,io.nsub); view(0,90);
[connectivity, coordinates, element_vals, ~] = makeParaviewData(xval, io.nsub, grafo);
data_struct = getSpectralData(io);
data_struct.force_type = parameters.force_type;
close all
save(strcat(filename,'.txt'),'connectivity','coordinates','element_vals','-ascii');
save(strcat(filename,'.mat'),'xval','x_history','fobj','fres', 'data_struct');
end