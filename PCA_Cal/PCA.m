function [mean_struct, Q_file_list, covariance, eigvec, eigval, total_no_frames] = PCA(trj_file_list, Q_trj_folder)
    [mean_struct, total_no_frames] = calc_mean_struct(trj_file_list);
    [Q_file_list] = write_Q_trj(trj_file_list, mean_struct, total_no_frames, Q_trj_folder);
    [covariance, eigvec, eigval] = calc_cov_n_decomposition(Q_file_list, numel(mean_struct));
end
