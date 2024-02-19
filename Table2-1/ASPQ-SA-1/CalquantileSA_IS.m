function [quantile_next] = CalquantileSA_IS(obser,p,stage,quantile_initial,ak,C,alpha)

%     step_size = ak / stage;
    step_size = ak / (C+ (stage)^(alpha));
    quantile_next = 0;
    n = size(obser,2);
    
    quantile_next = quantile_initial - step_size * ( (1 - (1/n) * ((obser(1,1:n,1) > quantile_initial)*obser(1,1:n,2)') ) - p );%%%%¥[¤Wobser LR      
end

