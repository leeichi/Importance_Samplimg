clear

p = [0.9];
k = [3];
n0 = [200 400 1000];
nb = [20 50 100];
trial = 500;
delta = 1;
mu = 0;
va = 100;
result = [];
c = 0.5;
v = 0.5;
count = 1;
for z = 1:length(k)
    for i = 1:length(n0)
        for j = 1:length(nb)
            [PCS,ANS ,CPU_TIME,quantile_TIME,dele_TIME,sample_TIME]=fsp_CMC_General(n0(i),p,k(z),c,v,nb(j),delta,mu,va,trial);%v=1/2
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