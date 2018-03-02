function sliceData = extractSlice(all_perm_data,col,L)
% for parallelization to eliminate large overhead transfering to and from workers
sliceData = all_perm_data.Data(col:L:end);