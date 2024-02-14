function  [quantile] = calQuantile(obs,s,L,n0,r,p,s_previously,amount)
    if(p>1)
        p = 1;
    end
    temp = squeeze(obs);
    temp_cdf = zeros(length(temp),1);
    
    for j = 1 : length(temp)%%%��C�Ӫ�cdf
        temp_cdf(j) = empCDF(temp(j),temp,s,r,L,n0,s_previously,amount);            
    end
    
    [~,idx] = sort(temp_cdf); % sort just the first column�Aidx���������Ȧb���Ӧ�m
    tmp_cdf = temp_cdf(idx);%%%%cdf�Ѥp�ƨ�j
    tmp = temp(idx,:);   % sort the whole matrix using the sort indices�A�̶��ǱƦC�Aquantile�ȥѤp�ƨ�j
    quantile = tmp(length(tmp) ,1);

    for j = 1 : length(tmp_cdf) 
        if(tmp_cdf(j) >= p)
            quantile = tmp(j);
            break
        end
    end
end

