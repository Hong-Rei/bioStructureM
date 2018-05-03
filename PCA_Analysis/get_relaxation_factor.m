function [relaxation_factor] = get_relaxation_factor(Normalized_AutoCorr,time_interval,correlation_time)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Get the relaxation factor (i.e. the diffusive part of time-correlation function)
%   The relaxation correlation function is modelled as an exponential function
% input:
%   Normalized_AutoCorr: Normalized Autocorrelation Function of a certain mode.
%   time_interval: The time elapse of each snapshot when you were saving MD trajectories. (in unit of ps)
%   diffusion_coefficient: The damping coefficient of Langevin dynamics for protein in solvent (in unit of 1/ps)
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	Num_of_frames = length(Normalized_AutoCorr);
	time_vector = time_interval*(0:(Num_of_frames-1));
	diffusion_coefficient = 1/correlation_time;
	exponent = ((-1)*diffusion_coefficient).*time_vector;
	relaxation_factor = exp(exponent);
	relaxation_factor = relaxation_factor./max(abs(relaxation_factor));

end

