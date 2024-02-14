function [psi_all_stage] = CalPsi(obser,p,stage,quantile_now,psi_all_previous,n0)


%------------------------------老師版本---------------------------------------    
    n = size(obser,2);
    psi_stage = (1 / n) * ((obser(1,1:n,1) > quantile_now) * (obser(1,1:n,2).^(2))') - ( 1 - p )^(2);
    
    if psi_stage < 0
        psi_stage = 0;    
    end
    
    if stage == 1 
        psi_all_stage = ( ((p*(1-p)) / n0) +  (psi_stage / n) );
    else
        psi_all_stage = (psi_all_previous + (psi_stage / n) );
    end
%------------------------------寬哥版本---------------------------------------
%     n = size(obser,2);
%     psi_stage = (1 / n) * ((obser(1,1:n,1) > quantile_initial) * (obser(1,1:n,2).^(2))') - ( 1 - p )^(2);
%     
%     if psi_stage < 0
%         psi_stage = 0;    
%     end
%     
%     if stage == 1
%         psi_all_stage = ((p*(1-p))/(n0) + psi_stage) / 2;    
%     else  
%     psi_all_stage = (stage / (stage + 1)) * psi_all_previous + (1 / (stage + 1)) * psi_stage;
%     end    
end

