function [max_PC_mode_to_exp_profile,correlation_information] = PCA_to_exp_profile_mapping_RMSD(exp_profile,reduced_profile)
%  Map the reduced PCA calculated profiles to the experimental profile using Root Mean-Squared Difference (RMSD).
% 
% Arguments:
%   exp_profile: Experimental profile.
%   reduced_profile: The (NxM) matrix containing the reduced profile calculated from PCA.
%   IMPORTANT: If the profile is NMR NOE order parameter, then order parameter profile should have missing residues removed!!
% 
% Returns:
%   max_PC_mode_to_exp_profile: The best mapped mode to the experimental profile. The matrix is (1x2).
%   The first column is the top mapped PC mode index, second is the RMSD.
%
%   correlation_information: The (Mx1) matrix containing mapping information.
%   The first column is the RMSD of each PC mode to experimental profile.
    num_of_modes = size(reduced_profile, 2);
    max_PC_mode_to_exp_profile = zeros(1,2);
    correlation_information = zeros(num_of_modes,1);

    for i = 1:num_of_modes
        current_rmsd = sqrt(mean((exp_profile - reduced_profile(:, i)).^2));
        correlation_information(i) = current_rmsd;
    end

    [min_RMSD,index] = min(correlation_information);
    max_PC_mode_to_exp_profile(1) = index;
    max_PC_mode_to_exp_profile(2) = min_RMSD;
end
