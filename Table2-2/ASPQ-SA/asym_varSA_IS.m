function[asym_variance] = asym_varSA_IS(p,b,w,r,quantile_initial_bw,psi_all_stage,stage)
    phi = zeros(1,length(b));
   
    for i = 1 : length(phi)
        if(p+b(i)>1)%%%%%%%%%%%%%%%?
            phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(29*(1-p)/15); %%初始階段如果很罕見 這邊會 = 0
        else    
            phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(2*b(i)); 
        end
    end


    Phi = w * phi';
    Psi = psi_all_stage * ( r / (stage)^(2) );
    asym_variance = (Psi*Phi.^2) / r;%%%%%有加.則不是內積

    
end

