function [PCA_to_GNM_mode_mapping] = PCA_to_GNM_mapping(reduced_Bval_from_PCA,reduced_Bval_from_GNM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	To do the frequency mapping from each GNM mode to PCA modes by comparing reduced B-factors.
%	Then the mapping results are saved.
% input:
%   reduced_Bval_from_PCA: The (NxM) matrix containing reduced B-factors from PCA.
%   reduced_Bval_from_GNM: The (NxM) matrix containing reduced B-factors from GNM.
% return:
%	PCA_to_GNM_mode_mapping: The mode mapping results. The matrix is (Mx4). M is the number of GNM modes.
%	The first column is the GNM mode index, second is the corresponding mapped PC mode index.
%	The third and fourth column contains the top correlation coefficient and the corresponding p-values (of the mapped PC mode). 
%
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[~,num_of_modes_PCA] = size(reduced_Bval_from_PCA);
	[~,num_of_modes_GNM] = size(reduced_Bval_from_GNM);
	PCA_to_GNM_mode_mapping = zeros(num_of_modes_GNM,4);
	PCA_to_GNM_mode_mapping(:,1) = (1:num_of_modes_GNM)';

	for i = 1:num_of_modes_GNM
		correlation_information = zeros(num_of_modes_PCA,2);
		for j = 1:num_of_modes_PCA
			[corr,p_value] = corrcoef(reduced_Bval_from_GNM(:,i),reduced_Bval_from_PCA(:,j));
			correlation_information(j,1) = corr(1,2);
			correlation_information(j,2) = p_value(1,2);
		end
		[max_corr,index] = max(correlation_information(:,1));
		PCA_to_GNM_mode_mapping(i,2) = index;
		PCA_to_GNM_mode_mapping(i,3) = max_corr;
		PCA_to_GNM_mode_mapping(i,4) = correlation_information(index,2);

		if PCA_to_GNM_mode_mapping(i,4) >= 0.05
			disp(['p_value is too high for GNM mode' num2str(i)]);
		end
	end

	save('PCA_to_GNM_mode_mapping.mat','PCA_to_GNM_mode_mapping','-v7.3');
end
