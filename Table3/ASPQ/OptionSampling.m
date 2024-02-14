function [payoff] = OptionSampling(S0,K,num,sigma,T,r,mu,var,t1)

dt=T;
T1=T-t1;
payoff = ones(num,1);%%%�h�֭Ӽ˥��x�}


for i = 1 : num %%%%�U�Ӫ�S
    
   

    a = normrnd(mu, var^(1/2));

    u=exp((r-1/2*sigma^2)*dt+sigma*sqrt(dt)*a);
    ut=exp((r-1/2*sigma^2)*t1+sigma*sqrt(t1)*a);

    S = S0*u;%%%%%% paper page.522������S����
    St = S0*ut;
    
    payoff(i) = blsprice(St,K,r,T1,sigma)-blsprice(S0,K,r,T,sigma);
end

end

