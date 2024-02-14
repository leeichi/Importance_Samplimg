clear

p = [0.99];
k = [10];
n0 = [600];
nb = [100];
trial = 500;
delta = 0.8;%[ak bk] = [300-1000]

mu = 0;
va = 1;%%%%§ï³o


result = [];
count = 1;


input_S0 = [105 102 98.9 95.75 92.65 89.45 86.25 83 79.7 76.4]; %%% initial price
input_K = [100 95 90 85 80 75 70 65 60 55];  %%% strike price

sigma = 0.2;
T = 6;
t1=0.5;
rate = 0.05;

for i = 1:length(n0)
    for j = 1:length(nb)
        for z = 1:length(k)
            [PCS,ANS ,CPU_TIME,quantile_TIME,dele_TIME,sample_TIME]=fsp_CMC_OptionSAA(n0(i),p,k(z),0.5,0.5,nb(j),delta,mu,va,trial,input_S0,input_K,sigma,T,rate,t1);
            result(count,1) = k(z);
            result(count,2) = n0(i);
            result(count,3) = nb(j);
            result(count,4) = PCS;
            result(count,5) = ANS;
            result(count,6) = CPU_TIME;
            result(count,7) = quantile_TIME;
            result(count,8) = dele_TIME;
            result(count,9) = sample_TIME;
            count = count + 1 ;
        end
    end
end