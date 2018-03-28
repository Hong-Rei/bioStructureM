%% get_ANM_period_power_law: Estimate the period of an ANM mode given its eigenvalue and the parameters of the power law
% 
% The power law has the following form: period = constant * eigenvalue.^exponent
% 
% Arguments:
% 	constant
% 	exponent
% 	eigenvalue
% 
% Returns:
% 	period (ns): Power law predicted period of ANM mode
function [period] = get_ANM_period_power_law(constant, exponent, eigenvalue)
	period = constant * eigenvalue.^exponent;
end