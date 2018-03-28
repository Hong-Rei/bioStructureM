function [RMSF] = getRMSFfromPCA(eigvalues,eigvectors,mode_selection)
%%%%%%%%%%%%%%%%%%%%need%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
%   eigvalues and eigvectors are from PCA covariance decomposition.
%   IMPORTANT: The input eigvalues and eigvectors should have smallest variance mode at top (MATLAB default!!) 
%               and the first six zero-modes should be already eliminated.
%   mode_selection is an array that specify the index of modes you want to use to reform covariance.
% return:
%   RMSF: an array which can be compared with experimental RMSF.
%   The unit of RMSF is in (angstorm).
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [three_N,num_of_modes] = size(eigvectors);
    num_of_atoms = three_N/3;

    if ~exist('mode_selection','var')
        mode_selection = 1:num_of_modes;
    end

    reduced_covariance = getReducedCovariancefromPCA(eigvalues,eigvectors,mode_selection);

    RMSF = sqrt(sum(reshape(diag(reduced_covariance), 3, num_of_atoms), 1));
end