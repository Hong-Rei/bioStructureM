%% get_ranks: Get the ranks of the values in ascending order.
% 
% Rank of identical values are assigned the average of their ranks
% 
% Argument:
%   vals(int/float): Vector of values to be ranked.
% 
% Returns:
%   ranks(float): Ranks of the values in ascending order.
function ranks=get_ranks(vals)
    unique_vals = unique(vals);
    [~,z1] = sort(vals);
    ranks = zeros(length(vals),1);
    ranks(z1) = 1:length(vals);

    for i=1:length(unique_vals)
        indices = unique_vals(i)==vals;
        ranks(indices,1) = mean(ranks(indices));
    end
end