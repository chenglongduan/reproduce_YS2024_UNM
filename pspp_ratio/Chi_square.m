function chi_misfit_norm = Chi_square(syn_vector, obs_vector, std_dev)

% Calculate Chi-squares misfit values


N_sample = length(syn_vector);
if N_sample ~= length(obs_vector)
    error('check number of samples in syn_vector and obs_vector');
end

chi_misfit = 0;

for i=1:N_sample

    chi_misfit = chi_misfit + ((syn_vector(i)-obs_vector(i))/std_dev(i))^2;

end

chi_misfit_norm = chi_misfit / N_sample;


end