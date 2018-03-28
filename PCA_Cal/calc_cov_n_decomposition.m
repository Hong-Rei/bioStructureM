%% Calculate the covariance matrix for PCA and decompose covariance matrix into it's eigenvector and eigenvalues
% 
%  Saves the covariance (covariance.mat), eigvec with non-zero eigval (eigvec.mat) and eigval (eigval.mat) to file
% 
%  Arguments:
%         no_of_trj_mat: Number of trajectories. File names starting with TRJ_new1.mat to TRJ_new(no_of_trj_mat).mat

function [covariance, eigvec, eigval] = calc_cov_n_decomposition(Q_file_list, N3)
    covariance = zeros(N3, N3);

    for traject_no = 1:numel(Q_file_list)
        Q_file = matfile(Q_file_list{traject_no});
        Q_file_variable_name = who(Q_file);
        Q_data = Q_file.(Q_file_variable_name{1});
        clear Q_file;

        covariance = covariance + Q_data * Q_data';
    end

    covariance = (covariance + covariance') / 2.0;

    [eigvec,eigval] = eig(covariance);
    eigval = diag(eigval);

    no_of_zeros = sum(eigval < 10e-6);
    disp(['No of zero eigvalues: ' num2str(no_of_zeros)])

    if no_of_zeros == 1
        eigvec = eigvec(:,2:end);
    else
        eigvec = eigvec(:,no_of_zeros + 1:end);
    end
end
