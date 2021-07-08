function outArray = chop_beats(beat_cellArray, t_vector)
%  Takes a cell array of beat-separated data and chops according to
%  shortest data vector; then, slaps on a time vector for easy keeping
%  track of
%  ------ Inputs ----------%
% beat_cellArray: Nx1 cellarray of variable length row vectors full of data
% t_vector: Nx1 time vector to store associated R peak times for each beat

% -------- Outputs --------%
% outArray: Nx(M+1) array of beat-separated data, where we have N beats and
% the chop length was M. We have an added column to keep track of
% time - this column is tacked on as the first column

% Initialize array to store beat lengths
length_beats = zeros(numel(beat_cellArray), 1);

% Store the lengths of each beat
for i = 1:length(length_beats)
    length_beats(i) = length(beat_cellArray{i});
end

% Use the shortest beat's length to chop
chop_length = min(length_beats);

% Initialize output array
outArray = zeros(length(length_beats), chop_length + 1);

% Store the time vector
outArray(:, 1) = t_vector;

% Store cell array data into output array, starting with the 2nd column
for i = 1:length(length_beats)
    for j = 2:chop_length + 1
        outArray(i, j) = beat_cellArray{i}(j-1);
    end
end

end

