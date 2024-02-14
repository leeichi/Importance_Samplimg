function[asym_variance] = asym_var(obs,b,s,L,n0,r,p,quantil,MU_Y,VAR_Y)
    obs = squeeze(obs);
    
    phi = zeros(1,length(b));
    for i = 1 : length(phi)
        if  p+b(i)>1 
            phi(i) = ((quantile(obs(:,1),1-(1-p)/30)) - quantile(obs(:,1),2*p-1+(1-p)/30))/(29*(1-p)/15);
        else
            phi(i) = (quantile(obs(:,1),p+b(i)) - quantile(obs(:,1),p-b(i)))/(2*b(i));

    
        end
    end
    Phi = phi';
    for i = 0: s-1
        if(i == 0)
            temp1 = obs((1:r),:);
            temp2(i+1) = ( 1 / length(temp1) ) * p * (1 - p);
        else
            temp1 = obs((n0+L*(i-1)+1:min(n0+L*i,r)),:);
            temp2(i+1) = ( ( 1 / ( length(temp1)^2 ) ) * ( (temp1(:,1)>quantil)' * temp1(:,2).^2 ) ) - ( ( 1 / length(temp1) ) * (1-p)^2);
        end
    
    end
   
    Psi = r * sum(temp2) / ((s+1)^2);%%看幾階除掉，意思為把前面的全部一起平均
    %%%%%%%%看要不要加/s^2，跟前面var一樣
    
    asym_variance = (Psi*Phi.^2)/r;

    clear temp1 temp2
    

    
end