function [all_names, indexMatrix]  = generateNameList(prefix, input_cell)
sz1 = length(input_cell);
idx_cell = cell(sz1,1);
msize = 1;
for i=1:sz1
    idx_cell{i} = 1:length(input_cell{i});
    msize = msize*length(idx_cell{i});
end

indices = cell(1, numel(idx_cell));   %output variable for ndgrid
[indices{:}] = ndgrid(idx_cell{:}); %get ngrid result

indexMatrix = zeros(msize,sz1);
for i=1:sz1
    array = indices{i};
    indexMatrix(:,i) = array(:);
end

all_names = cell(msize,1);
for i=1:msize
    midx = indexMatrix(i,:);
    str = prefix;
    for ii=1:sz1
        if ii == 1
            str = str + input_cell{ii}(midx(ii));
        else
            str = str +'_'+ input_cell{ii}(midx(ii));
        end
    end         
    all_names{i} = str;
end
all_names = string(all_names);
end