clear

delta = 0.8;%[ak bk] = [300-1000]
mu0 = 0;
va = 1;%%%%%100

trial =500;

p = [0.99];%%%%%%%%%%%%%%%調整
n0 = [600];
nb = [100];
k = [10];

ak = [800];%%%%%%%%%%%%%%%調整quantile
C = 0;
alpha = 1  ;

c = 0.5;
v = 0.5;

input_S0 = [105 102 98.9 95.75 92.65 89.45 86.25 83 79.7 76.4]; %%% initial price
input_K = [100 95 90 85 80 75 70 65 60 55];  %%% strike price



sigma = 0.2;
T = 6;
rate = 0.05;
t1=0.5;

result = [];
count = 1;

for z = 1:length(k)
    for i = 1:length(n0)
        for j = 1:length(nb)
                    [PCS,ANS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=ASPQ_IS_OptionSA(p,k(z),n0(i),nb(j),delta,mu0,va,trial,c,v,ak,C,alpha,input_S0,input_K,sigma,T,rate,t1);
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