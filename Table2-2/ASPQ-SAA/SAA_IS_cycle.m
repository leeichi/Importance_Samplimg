clear
mu = 0;
va = 100;
p = [0.9];
n0 = [200 400 1000];
nb = [20 50 100];
k = [3 10 50];
trial = 500;
count = 1;
result = [];
for z = 1:length(k)
    for i = 1:length(n0)
        for j = 1:length(nb)
               
            [ANS,PCS,CPU_TIME,quantile_TIME,optmize_TIME,dele_TIME,sample_TIME]=fsp_IS_quantile_batch(n0(i),p,k(z),nb(j),1,mu,va,trial);
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