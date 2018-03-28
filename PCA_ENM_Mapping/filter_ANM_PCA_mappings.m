%% filter_ANM_PCA_mappings: Filter ANM modes mapped to the same PCA mode. - Test
% 
% Multiple ANM modes can be mapped to the same PCA mode.
% The ANM mode with the highest correlation with a specific PCA mode will be 
% kept while the others will be filtered out.
% 
% Arguments:
%  ANM_PCA_mappings: Mapping between ANM and PCA where each row consist of [ANM mode, PCA mode, correlation]
% 
% Returns:
%  unique_ANM_PCA_mappings: Unique ANM-PCA mappings with the same data in each row as in ANM_PCA_mappings
function [unique_ANM_PCA_mappings] = filter_ANM_PCA_mappings(ANM_PCA_mappings)
    unique_PC_modes = unique(ANM_PCA_mappings(:,2));
    num_unique_PC_modes = size(unique_PC_modes, 1);
    unique_ANM_PCA_mappings = zeros(num_unique_PC_modes, size(ANM_PCA_mappings, 2));

    for idx = 1:num_unique_PC_modes
        PC_mode = unique_PC_modes(idx);
        select_ANM_modes = ANM_PCA_mappings(ANM_PCA_mappings(:,2) == PC_mode,:);

        [~, index_pos] = max(select_ANM_modes(:,3));
        unique_ANM_PCA_mappings(idx,:) = select_ANM_modes(index_pos,:);
    end
end