%% write_Q_trj: Align each frame to the mean structure and calculate Q
function [Q_file_list] = write_Q_trj(trj_file_list, mean_struct, total_no_frames, Q_trj_folder, mass)
    display('Writing Q files')

    no_of_trj_mat = numel(trj_file_list);
    Q_file_list = cell(no_of_trj_mat,1);

    trj_file = matfile(trj_file_list{1});
    trj_data_variable_name = who(trj_file);
    N3 = size(trj_file, trj_data_variable_name{1}, 1);
    N = N3/3;

    if ~exist('mass', 'var')
        mass = ones(N,3);
    end

    sum_of_mass = sum(mass(:,1));

    for traject_no = 1:no_of_trj_mat
        trj_file = matfile(trj_file_list{traject_no});
        trj_data_variable_name = who(trj_file);
        trj_data = trj_file.(trj_data_variable_name{1});
        clear trj_file;

        current_no_frames = size(trj_data,2);

        trj_data = reshape(trj_data, 3, N, current_no_frames);

        parfor frame_no = 1:current_no_frames
                 fromXYZ1 = bestfit(mean_struct, trj_data(:,:,frame_no)', mass, sum_of_mass);
                 trj_data(:,:,frame_no) = (fromXYZ1 - mean_struct)';
        end

        trj_data = reshape(trj_data, N3, current_no_frames) / sqrt(total_no_frames-1);
        save([Q_trj_folder 'Q_file_' num2str(traject_no) '.mat'], 'trj_data', '-v7.3');
        Q_file_list{traject_no} = [Q_trj_folder 'Q_file_' num2str(traject_no) '.mat'];

        clear trj_data;
    end
end