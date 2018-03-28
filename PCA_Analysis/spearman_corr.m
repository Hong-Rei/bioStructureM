%% spearman_corr: Calculate the Spearman correlation between 2 vectors of values
% 
% Arguments:
%   vals_1(int/float): Vector containing N values.
%   vals_2(int/float): Vector containing N values.
% 
% Returns:
%   corr_val(float): The Spearman correlation between vals_1 and vals_2
function [corr_val] = spearman_corr(vals_1, vals_2)
    rank_1 = get_ranks(vals_1);
    rank_2 = get_ranks(vals_2);

    [current_correl current_pval] = corrcoef(rank_1, rank_2);
    corr_val = current_correl(1,2);