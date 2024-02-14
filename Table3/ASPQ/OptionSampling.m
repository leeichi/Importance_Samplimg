function [payoff] = OptionSampling(S0,K,num,sigma,T,r,mu,var,t1)

dt=T;
T1=T-t1;
payoff = ones(num,1);%%%多少個樣本矩陣


for i = 1 : num %%%%各個的S
    
   

    a = normrnd(mu, var^(1/2));

    u=exp((r-1/2*sigma^2)*dt+sigma*sqrt(dt)*a);
    ut=exp((r-1/2*sigma^2)*t1+sigma*sqrt(t1)*a);

    S = S0*u;%%%%%% paper page.522中間的S公式
    St = S0*ut;
    
    payoff(i) = blsprice(St,K,r,T1,sigma)-blsprice(S0,K,r,T,sigma);
end

end

