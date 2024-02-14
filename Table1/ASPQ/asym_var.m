function[asym_variance] = asym_var(obs,b,r,p)
    obs = squeeze(obs);
    phi = zeros(1,length(b));
   
    for i = 1 : length(phi)
        if(p+b(i)>1)
            phi(i) = ((quantile(obs,1-(1-p)/30)) - quantile(obs,2*p-1+(1-p)/30))/(29*(1-p)/15); %%��l���q�p�G�ܨu�� �o��| = 0
        else    
            phi(i) = (quantile(obs,p+b(i)) - quantile(obs,p-b(i)))/(2*b(i)); 
        end
    end
    Phi = 1 * phi';
    Psi = p*(1-p);
    asym_variance = (Psi*Phi.^2)/r;%%%%%���[.�h���O���n
    
end