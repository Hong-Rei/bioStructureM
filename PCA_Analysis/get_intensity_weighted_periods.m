%% get_intensity_weighted_periods: Calculate the intensity weighted periods for all PCs
function [intensity_weighted_period] = get_intensity_weighted_periods(proj_file_list, traj_freq)
    num_of_proj_mat = numel(proj_file_list);
    proj_num_frames = {};
    proj_file_objs = {};
    proj_file_vars = {};

    for proj_num = 1:num_of_proj_mat
        proj_file_obj = matfile(proj_file_list{proj_num});
        proj_file_var_name = who(proj_file_obj);

        proj_file_objs{proj_num} = proj_file_obj;
        proj_num_frames{proj_num} = size(proj_file_obj, proj_file_var_name{1}, 2);
        proj_file_vars{proj_num} = proj_file_var_name{1};
    end

    total_num_frames = sum(cell2mat(proj_num_frames));
    proj_file_var_name = who(proj_file_objs{1});
    num_of_modes = size(proj_file_objs{1}, proj_file_var_name{1}, 1);
    num_of_modes_per_round = floor(num_of_modes/num_of_proj_mat);

    intensity_weighted_period = zeros(num_of_modes,1);
    temp_intensity_weighted_period = zeros(num_of_modes_per_round,1);
    current_traj = zeros(num_of_modes_per_round, total_num_frames, 'single');
    mode_end = 0;

    for iter_num = 1:(num_of_proj_mat-1)
        mode_start = mode_end + 1;
        mode_end = iter_num * num_of_modes_per_round;

        last_frame = 0;

        for proj_num = 1:num_of_proj_mat
            current_traj(:, last_frame + 1:last_frame + proj_num_frames{proj_num}) = proj_file_objs{proj_num}.(proj_file_vars{proj_num})(mode_start:mode_end, 1:proj_num_frames{proj_num});
            last_frame = last_frame + proj_num_frames{proj_num};
        end

        parfor mode_idx = 1:num_of_modes_per_round
            current_power_spectrum = double(get_PowerSpectrum(current_traj(mode_idx,:)));
            temp_intensity_weighted_period(mode_idx) = get_intensity_weighted_period(current_power_spectrum, traj_freq);
        end

        intensity_weighted_period(mode_start:mode_end) = temp_intensity_weighted_period;
    end

    clear current_traj;
    clear temp_intensity_weighted_period;

    remaining_modes = num_of_modes - mode_end;
    current_traj = zeros(remaining_modes, total_num_frames, 'single');
    temp_intensity_weighted_period = zeros(remaining_modes,1);
    mode_start = mode_end + 1;
    mode_end = num_of_modes;
    last_frame = 0;

    for proj_num = 1:num_of_proj_mat
        current_traj(:, last_frame + 1:last_frame + proj_num_frames{proj_num}) = proj_file_objs{proj_num}.(proj_file_vars{proj_num})(mode_start:mode_end, 1:proj_num_frames{proj_num});
        last_frame = last_frame + proj_num_frames{proj_num};
    end

    parfor mode_idx = 1:remaining_modes
        current_power_spectrum = double(get_PowerSpectrum(current_traj(mode_idx,:)));
        temp_intensity_weighted_period(mode_idx) = get_intensity_weighted_period(current_power_spectrum, traj_freq);
    end

    intensity_weighted_period(mode_start:mode_end) = temp_intensity_weighted_period;
end