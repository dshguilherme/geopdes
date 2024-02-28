function [missing_names, midx] = missingFileNames(str_prefix, all_names)
str_test = strcat(str_prefix,'*.mat');
files = dir(fullfile(pwd, str_test{1}));
fileNames = {files.name};
fileNames = strrep(fileNames, '.mat', '');

[missing_names, midx] = setdiff(all_names, fileNames);
end