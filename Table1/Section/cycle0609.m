clear
n0 = [200 400 1000];
p = [0.9];
k = [3 10 50];
nb = [20 50 100];
delta = 1;
va = 100;
trial = 500 ;
result = [];
count=1;
for i = 1:length(n0)
    for j = 1:length(nb)
        for z = 1:length(k)
            [PCS,ANS, CPU_TIME]=fsp_Sectioning(n0(i),p,k(z),nb(j),delta,va,trial);

            result(count,1) = k(z);
            result(count,2) = n0(i);
            result(count,3) = nb(j);
            result(count,4) = PCS;
            result(count,5) = ANS;
            result(count,6) = CPU_TIME;
          
            count = count + 1 ;
        end
    end
end

save('0609.mat')