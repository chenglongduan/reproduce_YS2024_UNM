function rms_value = RMSE(syn_vector, obs_vector)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


N_sample = length(syn_vector);
if N_sample ~= length(obs_vector)
    error('check number of samples in syn_vector and obs_vector');
end

res_sq = 0;

for i=1:N_sample
    res_sq = res_sq + (syn_vector(i)-obs_vector(i))^2;
end

rms_value = sqrt(res_sq / N_sample);


end