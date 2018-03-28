%% get_ANM_PCA_profile_mapping: Map ANM modes to PCA modes by profile matching
%  
% Each ANM mode is mapped to 1 PC mode by finding the profile with the highest
% Pearson correlation
% 
% Arguments:
%   reduced_ANM_profile: 3NxM real valued matrix where N are the number of residues
%                        and M is the m'th profile constructed using all ANM modes, m and larger. 
%   reduced_PCA_profile: 3NxK real valued matrix where N are the number of residues
%                        and K is the k'th profile constructed using all PCs, k and larger. 
% 
% Returns:
%   ANM_to_PCA_mapping: Mx4 real valued matrix where each m correspond to 1 ANM mode while the columns are
%                       ANM mode, PC mode, correlation and p-value.
function [ANM_to_PCA_mapping] = get_ANM_PCA_profile_mapping(reduced_ANM_profile, reduced_PCA_profile)
    total_num_ANM_modes = size(reduced_ANM_profile, 2);
    total_num_PC_mode = size(reduced_PCA_profile, 2);
    ANM_to_PCA_mapping = zeros(total_num_ANM_modes, 4);
    ANM_to_PCA_mapping(:, 1) = 1:total_num_ANM_modes;

    for ANM_idx = 1:total_num_ANM_modes
        ANM_profile = reduced_ANM_profile(:,ANM_idx);
        max_correl = 0.0;
        max_pval = 0.0;
        max_PC_idx = 0;
        
        for PC_idx = 1:total_num_PC_mode
            [correl_mat, pval] = corrcoef(ANM_profile, reduced_PCA_profile(:,PC_idx));
            current_correl = correl_mat(1,2);

            if current_correl > max_correl
                max_correl = current_correl;
                max_PC_idx = PC_idx;
                max_pval = pval(1,2);
            end
        end

        ANM_to_PCA_mapping(ANM_idx, 2) = max_PC_idx;
        ANM_to_PCA_mapping(ANM_idx, 3) = max_correl;
        ANM_to_PCA_mapping(ANM_idx, 4) = max_pval;
    end
end