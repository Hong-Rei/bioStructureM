%% est_ANM_IWP_power_law: Estimate the power law relationship between ANM eigenvalues and intensity weighted periods (IWP).
% 
% The power law: Est_Period (ns) = Constant * ANM_eigenvalue^(Exponent)
% 
% Arguments:
%  ANM_eigenvalues: Real valued vector containing M eigenvalues from M modes.
%  IWP: Contains M intensity weighted periods (ns) assigned to the respective ANM modes in ANM_eigenvalues.
% 
% Returns:
%  power_coefficient: Contains 2 values: [Exponent Constant]
%  correlation: The correlation between the IWP and power law predicted period
function [power_coefficient, correlation] = est_ANM_IWP_power_law(ANM_eigenvalues, IWP)
	power_coefficient = polyfit(log(ANM_eigenvalues), log(IWP), 1);
	power_coefficient(2) = exp(power_coefficient(2));

	correlation = corrcoef(get_ANM_period_power_law(power_coefficient(2), power_coefficient(1), ANM_eigenvalues), IWP);
	correlation = correlation(1,2);
end