addpath ../Autocorrelation_and_power_spectrum/
load auto_correlation_function.mat
[Num_of_modes,Num_of_frames] = size(auto_correlation_function);
characterstic_time = zeros(Num_of_modes,1);
relaxation_time = zeros(Num_of_modes,1);
for i = 1:Num_of_modes
    corr = auto_correlation_function(i,:);
    for j = 1:Num_of_frames
        if corr(j) <= 0.00001
            the_last_frame = j-1;
            break;
        end
    end
    current_corr = corr(1:the_last_frame);
    characterstic_time(i) = get_Characteristic_Time(current_corr,0.3);
    relaxation_time(i) = get_relaxation_time(current_corr,0.3);
end

save('characterstic_time.mat','characterstic_time','-v7.3');
save('relaxation_time.mat','relaxation_time','-v7.3');
