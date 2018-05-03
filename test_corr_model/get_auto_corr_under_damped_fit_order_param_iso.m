function [auto_corr_under_damped_fit_order_param_iso,B,rel_time,fluc_rel_time,fluc_time,resnorm,residual,exitflag,output] = get_auto_corr_under_damped_fit_order_param_iso(Normalized_AutoCorr,time_interval,init_B,init_rel_time,init_fluc_rel_time,init_fluc_time,options)


	if ~exist('init_B','var'),
    	init_B = 0.5;
	end

	if ~exist('init_rel_time','var')
    	init_rel_time = 5000;
	end

	if ~exist('init_fluc_rel_time','var')
    	init_fluc_rel_time = 5000;
	end

	if ~exist('init_fluc_time','var')
    	init_fluc_time = 10000;
	end


	Num_of_frames = length(Normalized_AutoCorr);
	time_vector = time_interval*(0:(Num_of_frames-1));

	lower_bound = [0,0,0,0];
	upper_bound = [1,100000,100000,1000000];
	init_params = [init_B,init_rel_time,init_fluc_rel_time,init_fluc_time,];
	fitting_auto_corr = @(params,time_vector)((exp((-1)*time_vector./params(2))).*((params(1)*cos(2*pi*(1/params(4))*time_vector))+((1-params(1)).*exp((-1)*time_vector./params(3)))));


	[fitted_params,resnorm,residual,exitflag,output] = lsqcurvefit(fitting_auto_corr,init_params,time_vector,Normalized_AutoCorr,lower_bound,upper_bound,options);
	B = fitted_params(1);
	rel_time = fitted_params(2);
	fluc_rel_time = fitted_params(3);
	fluc_time = fitted_params(4);
	auto_corr_under_damped_fit_order_param_iso = (exp((-1)*time_vector./rel_time)).*((B*cos(2*pi*(1/fluc_time)*time_vector))+((1-B).*exp((-1)*time_vector./fluc_rel_time)));

end
