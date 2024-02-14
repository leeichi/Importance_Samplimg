clear
mu = 0;
va = 1;

p = [0.99];

n0 = [600];
nb = [100];
k = [10];
trial = 500;
count = 1;
amount = 5;
result = [];
delta = 0.8 ;%[ak bk] = [300-1000]

input_S0 = [105 102 98.9 95.75 92.65 89.45 86.25 83 79.7 76.4]; %%% initial price
input_K = [100 95 90 85 80 75 70 65 60 55];  %%% strike price

sigma = 0.2;
T = 6;
rate = 0.05;
t1 = 0.5;

for i = 1:length(n0)
    for j = 1:length(nb)
        for z = 1:length(k)
               
            [ANS,PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_OptionSAA_LD(n0(i),p,k(z),nb(j),delta,mu,va,trial,input_S0,input_K,sigma,T,rate,t1,amount);
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
            count = count + 1 ;

        end
    end
end