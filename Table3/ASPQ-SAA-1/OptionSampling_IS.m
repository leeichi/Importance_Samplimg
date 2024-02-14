function [payoff,sample_path,LR_TIME] = OptionSampling_IS(S0,K,num,sigma,T,r,mu0,var0,mu,var,LR_TIME,t1)
dt=T;
T1=T-t1;

sample_path = [];
payoff = ones(num,1);%%%�h�֭Ӽ˥��x�}
for i = 1 : num %%%%�U�Ӫ�S
    
    LR = 1;
    a = normrnd(mu, var^(1/2));
    
    timer6=tic;
    LR = LR * (normpdf(a,mu0,var0^(1/2)) ./ normpdf(a,mu,var^(1/2)));
    LR_TIME = LR_TIME + toc(timer6);

    u=exp((r-1/2*sigma^2)*dt+sigma*sqrt(dt)*a);
    ut=exp((r-1/2*sigma^2)*t1+sigma*sqrt(t1)*a);

    S = S0*u;%%%%%% paper page.522������S����
    St = S0*ut;

    
    payoff(i,1) = blsprice(St,K,r,T1,sigma)-blsprice(S0,K,r,T,sigma);
    timer6=tic;
    payoff(i,2) = LR;
    sample_path = [sample_path;a];
    LR_TIME = LR_TIME + toc(timer6);
    
end

end

