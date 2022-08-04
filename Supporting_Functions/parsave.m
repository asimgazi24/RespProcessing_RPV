function parsave(fname, data_struct)
% This function takes as input a file name and a struct with fields
% containing all the data that needs to be saved in a parfor loop. This
% function saves each field on its own using its field name

% ------------ Inputs ---------------- %
% fname: string storing the desired filename to save to
% data_struct: struct with fieldnames corresponding to what the variables
% should be named when saving the data to the file

save(fname, '-struct', 'data_struct');

end

