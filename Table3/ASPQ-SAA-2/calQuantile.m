function  [quantile] = calQuantile(obs,s,L,n0,r,p,s_previously,amount)
    if(p>1)
        p = 1;
    end
    temp = squeeze(obs);
    temp_cdf = zeros(length(temp),1);
    
    for j = 1 : length(temp)%%%算每個的cdf
        temp_cdf(j) = empCDF(temp(j),temp,s,r,L,n0,s_previously,amount);            
    end
    
    [~,idx] = sort(temp_cdf); % sort just the first column，idx為對應的值在哪個位置
    tmp_cdf = temp_cdf(idx);%%%%cdf由小排到大
    tmp = temp(idx,:);   % sort the whole matrix using the sort indices，依順序排列，quantile值由小排到大
    quantile = tmp(length(tmp) ,1);

    for j = 1 : length(tmp_cdf) 
        if(tmp_cdf(j) >= p)
            quantile = tmp(j);
            break
        end
    end
end

