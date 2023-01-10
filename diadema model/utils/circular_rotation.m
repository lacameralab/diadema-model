function [rotated_data] = circular_rotation(data,shift)
% [rotated_data] = circular_rotation(data,shift)
% Circular rotate the matrix data in row. Shift the matrix with the number
% of column "shift". shift = 2, rotate the matrix by 2. The first column
% moves to the third column. The last column moves to the second column.
% In general, the ith column moves to the 1+mod(i+shift,ncol) column. 
% Input:
%   data: spike count matrix (e.g. spkc{1}) for pHMM, sequence matrix (e.g.
%         seq{1}) for bHMM.
% 
% Tianshu Li, 2021

ncol = size(data,2);

new_indx = 1+mod((1:ncol)-shift-1,ncol);
rotated_data = data(:,new_indx);

end