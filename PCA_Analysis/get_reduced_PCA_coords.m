%% get_reduced_PCA_coords: Generate the reduced PCA coords from projection vector
% 
% The coordinate differences between a structure and the PCA's superposition structure
% (after remove translational & rotational differences) can be rebuilt from the projection
% vector.
% 
% This is calculated by scaling each eigenvector by the magnitude of the projection and then
% summing up these scaled eigenvectors.
% 
% The projection vector contains the projection of the coordinate differences to each PCA's
% eigenvector.
% 
% The reduced PCA coords are coordinate differences rebuilt after removing the first k-1 PCs.
% 
% Arguments:
%   eig_vects (MxN single/double matrix): The N columns represent the eigenvectors
%       while M = 3 x num of atoms. The eigenvectors should be ordered such that PC1
%       is located on the RIGHT most column while the last PC is on the LEFT most column.
% 
%   proj_vect (Nx1 single/double vector): The N projections of the coordinate differences to
%       the N PC's. Projections should be ordered starting from the last PC's projection to PC1's
%       projection.
% 
% Returns:
%   reduced_coords (MxN single/double matrix): The N columns of the matrix represent the reduced
%   	PCA coordinates starting with k=1 to k=number of eigenvectors.
function [reduced_coords] = get_reduced_PCA_coords(eig_vects, proj_vect)
    scaled_eig_vects = bsxfun(@times, eig_vects, proj_vect');
    reduced_coords = cumsum(scaled_eig_vects, 2);
end