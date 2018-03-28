function [mean_struct, total_no_frames] = calc_mean_struct(trj_file_list, mass)
    no_of_trj_mat = numel(trj_file_list);

    % Load first trajectory and use first frame as reference structure
    trj_file = matfile(trj_file_list{1});
    trj_data_variable_name = who(trj_file);

    N3 = size(trj_file, trj_data_variable_name{1}, 1);
    N = N3/3;
    mean_struct = double(reshape(trj_file.(trj_data_variable_name{1})(:,1),3,N)');
    clear trj_file;

    if ~exist('mass', 'var')
        mass = ones(N,3);
    end

    % Move COM to 0,0,0
    sum_of_mass = sum(mass(:,1));
    mass_center = sum(mean_struct .* mass) ./ sum_of_mass;
    mean_struct = mean_struct - repmat(mass_center,N,1);

    checkrmsd = 1e10;
    total_no_frames = 0;

    while checkrmsd >= 1e-6
        current_structure = zeros(size(mean_struct));
        total_no_frames = 0;

        % Loop through all trajectories and sum the aligned structures to current_structure 
        for traject_no = 1:no_of_trj_mat
            trj_file = matfile(trj_file_list{traject_no});
            trj_data_variable_name = who(trj_file);

            current_no_frames = size(trj_file, trj_data_variable_name{1}, 2);
            total_no_frames = total_no_frames + current_no_frames;

            trj_data = reshape(trj_file.(trj_data_variable_name{1}), 3, N, current_no_frames);
            clear trj_file;

            parfor frame_no = 1:current_no_frames
                fromXYZ1 = bestfit(mean_struct, trj_data(:,:,frame_no)', mass, sum_of_mass);

                % Make sure current_structure is a double float since we're summing up the structures
                current_structure = current_structure + double(fromXYZ1);
            end

            clear trj_data;
        end

        [~, ~, checkrmsd, ~, mean_struct] = rmsdfitm(mean_struct, current_structure/total_no_frames, mass(:,1));
        display(['Current Round:' num2str(checkrmsd)])
    end
end
