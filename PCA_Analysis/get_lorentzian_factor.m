function [lorentzian_factor] = get_lorentzian_factor(Normalized_PowerSpectrum,correlation_time,time_interval)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Get the relaxation factor (i.e. the diffusive part of time-correlation function)
%   The relaxation correlation function is modelled as an exponential function
% input:
%   Normalized_AutoCorr: Normalized Autocorrelation Function of a certain mode.
%   time_interval: The time elapse of each snapshot when you were saving MD trajectories. (in unit of ps)
%   diffusion_coefficient: The damping coefficient of Langevin dynamics for protein in solvent (in unit of 1/ps)
% Editor: Hong Rui
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	half_Num_of_frames = length(Normalized_PowerSpectrum);
	lorentzian_factor = zeros(1,half_Num_of_frames);
	sampling_frequency = 1/time_interval;
	frequency_vector = sampling_frequency*(0:(half_Num_of_frames-1))./(half_Num_of_frames*2);
	diffusion_coefficient = 1/correlation_time;

	for i = 1:half_Num_of_frames
		lorentzian_factor(i) = (10^12)*diffusion_coefficient/((diffusion_coefficient^2)+(4*pi^2*frequency_vector(i)^2));
	end
	
end