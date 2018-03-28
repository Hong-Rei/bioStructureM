%% project_frames_to_PCs: Project each frame to each PC
% 
%  Arguments:
%         eigvec: Matrix of (No of atoms x 3) by (No of PCs)
%         frames: Frames of coordinate of atoms aligned to the mean structure used in PCA.
%                 (No of atoms * 3) by (No of frames) matrix
%  Returns:
%         projections: Matrix of (No of PCs) by (No of frames)

function [proj_file_list] = project_frames_to_PCs(Q_file_list, eigvec, total_no_frames)

    % Assuming eigenvectors are ordered based on ascending order of eigenvalues
    % Switch it to descending order and then transpose eigvec
    eigvec = eigvec(:,end:-1:1)';

    no_Q_files = numel(Q_file_list);
    constant = sqrt(total_no_frames - 1);
    proj_file_list = cell(no_Q_files,1);

    for Q_file_no = 1:no_Q_files
        Q_file = matfile(Q_file_list{Q_file_no});
        Q_file_variable_name = who(Q_file);

        current_projection = eigvec * Q_file.(Q_file_variable_name{1}) * constant;
        save(['proj_' num2str(Q_file_no) '.mat'], 'current_projection', '-v7.3');
        proj_file_list{Q_file_no} = ['proj_' num2str(Q_file_no) '.mat'];
    end
end
