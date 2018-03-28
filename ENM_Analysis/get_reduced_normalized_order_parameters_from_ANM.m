function [reduced_normalized_order_parameter_from_ANM] = get_reduced_normalized_order_parameters_from_ANM(PDB_Structure,eigvalues)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To calculate the normalized reduced NMR NOE order parameters (S^2) by stepwise removing the slowest mode contribution.
% 
% The normalization ensures that the sum of the variances along the x, y and z dimensions of the bond equals to 1
% 
% input:
%   eigvalues is the ANM eigenvalues after removing six zero modes.
%   IMPORTANT!! The eigvalues and PDB_Structure should only contain NH atom ONLY!!
% return:
%   reduced_normalized_order_parameter_from_ANM: A (NXM) matrix, N is the residue number and M is the mode number. 
%
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eigvectors = getANM(PDB_Structure);
    [SIX_N,num_of_modes] = size(eigvectors);
    num_of_atoms = SIX_N/3;
    num_of_res = num_of_atoms/2;    
    reduced_normalized_order_parameter_from_ANM = zeros(num_of_res,num_of_modes);

    for i = 1:num_of_modes
        mode_selection = i:num_of_modes;
        reduced_normalized_order_parameter_from_ANM(:,i) = getNormalizedNMROrderParameterfromANM(PDB_Structure,eigvalues,mode_selection);
    end

end
