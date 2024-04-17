function std_1col = amp_ratio_std_bsp(numer_datamat,denom_datamat,geom_tab_new,dt,t_numer,t_denom,twin_numer_prior,twin_numer_post,...
    twin_denom_prior,twin_denom_post,bin_min,bin_incre,bin_max,method,bsp_Nrepeat)

%   Bootstrap resampling to calculate standard deviation (i.e. uncertainty)
%   Resample all the traces and repeat this process for N times
%   bsp_Nrepeat means bootstrap repeating times (>=500, 1000 is good)

ntr_new = length(t_numer);
if ntr_new ~= size(numer_datamat,2)
    error('wrong trace number!');
end

bins = bin_min:bin_incre:bin_max;
nbins = length(bins);
ratio_store = zeros(nbins-1,bsp_Nrepeat);
std_1col = zeros(nbins-1,1);

for ibsp=1:bsp_Nrepeat
    
    trace_no_tmp = randi([1 ntr_new],ntr_new,1); % random number generator
    trace_no = sort(trace_no_tmp);

    numer_data_bsp = numer_datamat(:,trace_no);
    denom_data_bsp = denom_datamat(:,trace_no);
    geom_tab_bsp = geom_tab_new(trace_no,:);
    t_numer_bsp = t_numer(trace_no);
    t_denom_bsp = t_denom(trace_no);

    % main body
    amp_numer_tmp = zeros(ntr_new,1);
    amp_denom_tmp = zeros(ntr_new,1);
    ratio_tmp = zeros(ntr_new,1);

    for i=1:ntr_new
        it1_numer = floor((t_numer_bsp(i)-twin_numer_prior)/dt);
        it2_numer = ceil((t_numer_bsp(i)+twin_numer_post)/dt);
    
        it1_denom = floor((t_denom_bsp(i)-twin_denom_prior)/dt);
        it2_denom = ceil((t_denom_bsp(i)+twin_denom_post)/dt);
        
        if strcmp(method,'max')
            amp_numer_tmp(i) = max(numer_data_bsp(it1_numer:it2_numer,i));
            amp_denom_tmp(i) = max(denom_data_bsp(it1_denom:it2_denom,i));
        end
        if strcmp(method,'median')
            amp_numer_tmp(i) = median(numer_data_bsp(it1_numer:it2_numer,i));
            amp_denom_tmp(i) = median(denom_data_bsp(it1_denom:it2_denom,i));
        end
        if strcmp(method,'mean')
            amp_numer_tmp(i) = mean(numer_data_bsp(it1_numer:it2_numer,i));
            amp_denom_tmp(i) = mean(denom_data_bsp(it1_denom:it2_denom,i));
        end
    
        ratio_tmp(i) = amp_numer_tmp(i)/amp_denom_tmp(i);
    end

    for j=1:nbins-1
        if j==1
            index = find(geom_tab_bsp(:,5)>=bins(j) & geom_tab_bsp(:,5)<=bins(j+1));
        else
            index = find(geom_tab_bsp(:,5)>bins(j) & geom_tab_bsp(:,5)<=bins(j+1));
        end
    
        ratio_store(j,ibsp) = median(ratio_tmp(index)); % mean, median
    end

end

% standard deviation
for ibin=1:nbins-1
    std_1col(ibin) = std(ratio_store(ibin,:),1);
end


end