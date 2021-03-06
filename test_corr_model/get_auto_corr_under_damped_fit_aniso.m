function [auto_corr_under_damped_fit_aniso,A,rel_time_1,rel_time_2,fluc_time,resnorm,residual,exitflag,output] = get_auto_corr_under_damped_fit_aniso(Normalized_AutoCorr,time_interval,init_A,init_rel_time_1,init_rel_time_2,init_fluc_time,options)

	if ~exist('init_A','var')
    	init_A = 0.5;
	end

	if ~exist('init_rel_time_1','var')
    	init_rel_time_1 = 5;
	end

	if ~exist('init_rel_time_2','var')
    	init_rel_time_2 = 5;
	end

	if ~exist('init_fluc_time','var')
    	init_fluc_time = 10;
	end


	Num_of_frames = length(Normalized_AutoCorr);
	time_vector = time_interval*(0:(Num_of_frames-1));

	lower_bound = [0,0,0,0];
	upper_bound = [1,100000,100000,1000000];
	init_params = [init_A,init_rel_time_1,init_rel_time_2,init_fluc_time];
	fitting_auto_corr = @(params,time_vector)(((params(1)*exp((-1)*time_vector./params(2)))+((1-params(1))*exp((-1)*time_vector./params(3)))).*cos(2*pi*(1/params(4))*time_vector));


	[fitted_params,resnorm,residual,exitflag,output] = lsqcurvefit(fitting_auto_corr,init_params,time_vector,Normalized_AutoCorr,lower_bound,upper_bound,options);
	A = fitted_params(1);
	rel_time_1 = fitted_params(2);
	rel_time_2 = fitted_params(3);
	fluc_time = fitted_params(4);
	auto_corr_under_damped_fit_aniso = ((A*exp((-1)*time_vector./rel_time_1))+((1-A)*exp((-1)*time_vector./rel_time_2))).*cos(2*pi*(1/fluc_time)*time_vector);

end