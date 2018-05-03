function [Normalized_AutoCorr_remove_rel] = get_Normalized_AutoCorr_remove_rel(Normalized_AutoCorr,relaxation_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Get the relaxation factor (i.e. the diffusive part of time-correlation function)
%   The relaxation correlation function is modelled as an exponential function
% input:
%   Normalized_AutoCorr: Normalized Autocorrelation Function of a certain mode.
%   time_interval: The time elapse of each snapshot when you were saving MD trajectories. (in unit of ps)
%   diffusion_coefficient: The damping coefficient of Langevin dynamics for protein in solvent (in unit of 1/ps)
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Normalized_AutoCorr_remove_rel = Normalized_AutoCorr - relaxation_factor;
	Normalized_AutoCorr_remove_rel = Normalized_AutoCorr_remove_rel./max(abs(Normalized_AutoCorr_remove_rel));
	
end
