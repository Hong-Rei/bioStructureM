%% calc_reduced_order_param: Calculate the reduced order parameter from the projected trajectory.
% 
% The trajectory contains the projection of each frame's coordinate differences onto each PCA mode.
% The coordinate differences are between each frame's coordinate and the PCA's superposition structure
% (after remove translational & rotational differences).
% 
% Arguments:
%   trj_file_list (Nx1 cell array): Contains the full file path to the N projection trajectory files.
% 
%   N_index (Nx1 int32/int64 array): The index (starts at 1) position of the N-backbone atoms.
% 
%   H_index (Nx1 int32/int64 array): The index (starts at 1) position of the hydrogens bonded to the N-backbone atoms.
% 
%   eig_vects (MxN single/double matrix): The N columns represent the eigenvectors while M = 3 x num of atoms.
%       The eigenvectors should be ordered such that PC1 is located on the LEFT most column while the last PC is on
%       the RIGHT most column.
% 
%   mean_struct (3Nx1 single/double vector): The mean/superposition structure's coordinates containing N atoms.
%       The coordinates are ordered: x1, y1, z1, x2, y2, z2, ..., xN, yN, zN.
% 
% Returns:
%   reduced_order_param (NxM double vector): The order parameters for N residues after removing the first k-1
%       eigenvectors, where k=1 to k=M. The left most column is the order parameter calculated using all M PCs
%       while the right most column's order parameter was calculated using the last PC.
function [reduced_order_param] = calc_reduced_order_param(trj_file_list, N_index, H_index, eigvects, mean_struct)
    no_of_trj_mat = numel(trj_file_list);
    total_num_frames = 0;
    eigvects = fliplr(eigvects);

    num_PC_modes = size(eigvects, 2);
    num_residues = numel(N_index);

    sqrd_NH_vect_traj = zeros(num_residues * 3, num_PC_modes, 'double');
    xy_NH_vect_traj = zeros(num_residues, num_PC_modes, 'double');
    xz_NH_vect_traj = zeros(num_residues, num_PC_modes, 'double');
    yz_NH_vect_traj = zeros(num_residues, num_PC_modes, 'double');

    N_3N_start_index = (N_index-1) * 3;
    N_atom_select = [(N_3N_start_index+1) (N_3N_start_index+2) (N_3N_start_index+3)]';
    N_sel_eigvects = double(eigvects(N_atom_select,:));

    H_3N_start_index = (H_index-1) * 3;
    H_atom_select = [(H_3N_start_index+1) (H_3N_start_index+2) (H_3N_start_index+3)]';
    H_sel_eigvects = double(eigvects(H_atom_select,:));

    NH_mean_diff = double(mean_struct(H_atom_select,:) - mean_struct(N_atom_select,:));

    for traj_idx = 1:no_of_trj_mat
        trj_file = matfile(trj_file_list{traj_idx});
        trj_data_variable_name = who(trj_file);

        current_num_frames = size(trj_file, trj_data_variable_name{1}, 2);
        total_num_frames = total_num_frames + current_num_frames;
        current_proj_trj = flipud(trj_file.(trj_data_variable_name{1})(1:num_PC_modes,1:current_num_frames));

        parfor current_frame = 1:current_num_frames
            current_proj_vect = double(current_proj_trj(:,current_frame));

            N_reduced_coords = get_reduced_PCA_coords(N_sel_eigvects, current_proj_vect);
            H_reduced_coords = get_reduced_PCA_coords(H_sel_eigvects, current_proj_vect);

            NH_reduced_vects = bsxfun(@plus, (H_reduced_coords - N_reduced_coords), NH_mean_diff); % Vector from N to H
            vects_magnitude = sqrt(sum(reshape(NH_reduced_vects, 3, num_residues, num_PC_modes).^2, 1));
            NH_reduced_vects_normed = reshape(bsxfun(@rdivide, reshape(NH_reduced_vects, 3, num_residues, num_PC_modes), vects_magnitude), num_residues * 3, num_PC_modes);

            sqrd_NH_vect_traj =  sqrd_NH_vect_traj + NH_reduced_vects_normed.^2;
            xy_NH_vect_traj = xy_NH_vect_traj + NH_reduced_vects_normed(1:3:end, 1:end) .* NH_reduced_vects_normed(2:3:end, 1:end);
            xz_NH_vect_traj = xz_NH_vect_traj + NH_reduced_vects_normed(1:3:end, 1:end) .* NH_reduced_vects_normed(3:3:end, 1:end);
            yz_NH_vect_traj = yz_NH_vect_traj + NH_reduced_vects_normed(2:3:end, 1:end) .* NH_reduced_vects_normed(3:3:end, 1:end);
        end
    end

    reduced_order_param = 3/2 * (squeeze(sum(reshape((sqrd_NH_vect_traj./total_num_frames).^2, 3, num_residues, num_PC_modes), 1)) + 2*((xy_NH_vect_traj./total_num_frames).^2 + (xz_NH_vect_traj./total_num_frames).^2 + (yz_NH_vect_traj./total_num_frames).^2)) - 1/2;
    reduced_order_param = fliplr(reduced_order_param);
end