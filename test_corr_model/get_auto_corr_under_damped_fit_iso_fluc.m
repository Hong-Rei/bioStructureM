function [auto_corr_under_damped_fit_iso_fluc,fluc_time,resnorm,residual,exitflag,output] = get_auto_corr_under_damped_fit_iso_fluc(Normalized_AutoCorr,time_interval,rel_time,init_fluc_time,options)

	if ~exist('init_fluc_time','var')
    	init_fluc_time = 10;
	end

	Num_of_frames = length(Normalized_AutoCorr);
	time_vector = time_interval*(0:(Num_of_frames-1));

	lower_bound = 0;
	upper_bound = 1000000;
	init_time_scale = init_fluc_time;
	fitting_auto_corr = @(time_scale,time_vector)(exp((-1)*time_vector./rel_time).*cos(2*pi*(1/time_scale).*time_vector));


	[time_scale,resnorm,residual,exitflag,output] = lsqcurvefit(fitting_auto_corr,init_time_scale,time_vector,Normalized_AutoCorr,lower_bound,upper_bound,options);
	fluc_time = time_scale;
	auto_corr_under_damped_fit_iso_fluc = exp((-1)*time_vector./rel_time).*cos(2*pi*(1/fluc_time).*time_vector);

end