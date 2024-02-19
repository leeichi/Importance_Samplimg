clear

delta = 1; % indifferent zone parameter
mu0 = 0; % mean
va = 100; % vatiance

trial = 500;  %% trial

p = [0.99]; % p-quantile
n0 = [200 400 1000]; % initial stage sample size
nb = [20 50 100]; % other stage sample size
k = [3 10 50]; % number of system
ak = [800]; % step-size for update quantile

% step-size parameter
C = 0;
alpha = 1;

% bandwidth parameter
c = 0.5;
v = 0.5;
% norminv(p,mu0,va^(1/2))
result = [];
count = 1;

for z = 1:length(k)
    for i = 1:length(n0)
        for j = 1:length(nb)
                    [PCS,ANS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=ASPQ_SA_IS(p,k(z),n0(i),nb(j),delta,mu0,va,trial,c,v,ak,C,alpha);
                    result(count,1) = k(z);
                    result(count,2) = n0(i);
                    result(count,3) = nb(j);
                    result(count,4) = PCS;
                    result(count,5) = ANS;
                    result(count,6) = CPU_TIME;
                    result(count,7) = quantile_TIME;
                    result(count,8) = optmize_TIME;
                    result(count,9) = dele_TIME;
                    result(count,10) = sample_TIME;
                    
                    count = count + 1;
        end
    end
end