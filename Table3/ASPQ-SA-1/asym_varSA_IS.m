function[asym_variance] = asym_varSA_IS(p,b,w,r,quantile_initial_bw,psi_all_stage,stage)
    phi = zeros(1,length(b));
   
    for i = 1 : length(phi)
        if(p+b(i)>1)%%%%%%%%%%%%%%%?
            phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(29*(1-p)/15); %%初始階段如果很罕見 這邊會 = 0
            %phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(9*(1-p)/5); %%初始階段如果很罕見 這邊會 = 0
%            phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(2*b(i));
%             phi(i) = (CalquantileSA(obser,1-(1-p)/10,stage,quantile_last,a) - CalquantileSA(obser,2*p-1+(1-p)/10,stage,quantile_last,a))/(9*(1-p)/5);
        else    
            phi(i) = (quantile_initial_bw(1) - quantile_initial_bw(2))/(2*b(i)); 
%             phi(i) = (CalquantileSA(obser,p+b(i),stage,quantile_last,a) - CalquantileSA(obser,p-b(i),stage,quantile_last,a))/(2*b(i));
        end
    end
    Phi = w * phi';
% CalquantileSA(obser,p,stage,quantile_last,a);
%% 
%      Psi = p*(1-p);
     Psi = psi_all_stage * ( r / (stage)^(2) );
%      Psi = ( ((obser(1,1:n,1) > quantile_initial) * (obser(1,1:n,2).^(2))' / n) - ( 1 - p )^(2));
%      psi_all_stage
     
     asym_variance = (Psi*Phi.^2) / r;%%%%%有加.則不是內積
%      asym_variance = (Psi*Phi.^2) / (r * (stage^(2)));%%%%%有加.則不是內積
    
end

